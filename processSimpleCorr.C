/* 
 * This macro compares our results to published
 * http://hepdata.cedar.ac.uk/view/ins930312
 * without correction in vertex bins, without mixed event correction
 * with TPCOnly cut, in pTt 8-15, pTa 3-4, 4-6, 6-8, 8-10 GeV
 */

// TODO add guards
#include "../inc/AliJHistManagerROOT6.cxx"
#include "../inc/mfit.h"
#include "../inc/mtools.h"

// -------------------------------------
// Simplified correlation analysis class
//
class MSimpleCorr
{
    private:
        TFile * fInFile;

    public:
        AliJBin *fCent, *fVtx, *fPTt, *fPTa, *fEta, *fPhi;  
        int fNumPhi, fNumEta, fNumPta, fNumPtt, fNumVtx, fNumCent;
        int fVertexSkip, fPhiSkip;

        AliJHistManager * fhst;
        AliJTH1D fhDEtaNearRaw, fhDEtaNearMix, fhTriggPt, fhZVert;

        double fPtaPublishedBins[6];
        const int fNumPtaPublished = 5;

        TString type = "";
        TString TrackCut = "";
        TString fOutFileName;

        // INITIALIZATION
        MSimpleCorr(TString inName, TString trackC, double vskip, double phiskip) 
        {
            fPtaPublishedBins[0]=2;fPtaPublishedBins[1]=3;fPtaPublishedBins[2]=4;fPtaPublishedBins[3]=6;fPtaPublishedBins[4]=8;fPtaPublishedBins[5]=10;
            TrackCut = trackC;

            fInFile = TFile::Open(inName);
            fInFile->cd(); //InFile->ls();

            // hist manager and histos from input file
            fhst = new AliJHistManager( "hst" );
            fhst->LoadConfig( Form("JDiHadronCorr_%s/AliJHistos", TrackCut.Data()) );
            fhDEtaNearRaw = fhst->GetTH1D("hDEtaNear");           // 5D: C, V, E, T, A
            fhDEtaNearMix = fhst->GetTH1D("hDEtaNearM");          // 5D: C, V, E, T, A
            fhTriggPt     = fhst->GetTH1D("hTriggPtBin");         // 3D: C, V, T
            fhZVert       = fhst->GetTH1D("hZVert");              // 1D: C

            // binning 
            fCent = fhst->GetBin("Cent");   fNumCent = fCent->Size();
            fVtx  = fhst->GetBin("Vtx");    fNumVtx  = fVtx->Size();
            fPTt  = fhst->GetBin("PTt");    fNumPtt  = fPTt->Size();
            fPTa  = fhst->GetBin("PTa");    fNumPta  = fPTa->Size();
            fEta  = fhst->GetBin("EtaGap"); fNumEta  = fEta->Size();
            fPhi  = fhst->GetBin("PhiGap"); fNumPhi  = fPhi->Size();

            fVertexSkip = fVtx->GetBin( vskip );
            fPhiSkip    = fPhi->GetBin( phiskip ); // it is always near-side, so 1.6 can be used here

            if( inName.Contains("pp") )   type = "pp";
            if( inName.Contains("PbPb") ) type = "PbPb";

            type=="PbPb" ? fNumCent = 2 : fNumCent=1;

            cout << "initialization done..." << endl;
            LoadNtrigg();
        }

        // DESTRUCTOR
        ~MSimpleCorr() 
        {
            for(int iptt=1; iptt<fNumPtt; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if( fPTt->At(iptt) < fPTa->At(ipta) )
                            continue;
                        delete hDEtaRaw[ic][iptt][ipta], hDEtaMix[ic][iptt][ipta], hDEtaReal[ic][iptt][ipta];
                    }
                    for(int ipta_new=0; ipta_new<fNumPtaPublished; ipta_new++){
                        if( fPTt->At(iptt) < fPtaPublishedBins[ipta_new] )
                            continue;
                        delete hDEtaRawNew[ic][iptt][ipta_new], hDEtaRealNew[ic][iptt][ipta_new], hDEtaMixNew[ic][iptt][ipta_new], hDEtaRawFlip[ic][iptt][ipta_new], hDEtaRealFlip[ic][iptt][ipta_new];
                    }
                    for(int ifit=0; ifit<kF; ifit++){
                        delete hYield[ifit][ic][iptt], hYield_Int[ifit][ic][iptt];
                    }
                }
            }
            fInFile->Close();
        }



        const static int kT=5, kA=10, kC=5, kS=3, kF=3;
        // correlation histos: eta 
        TH1D * hDEtaRaw[kC][kT][kA], * hDEtaMix[kC][kT][kA], * hDEtaReal[kC][kT][kA];

        // output histos with new pta binning
        TH1D * hYield[kF][kC][kT];
        // yield calculated by integrating histograms. Constans estimated from the different fits.
        TH1D * hYield_Int[kF][kC][kT];  

        TH1D * hDEtaRawNew[kC][kT][kA], * hDEtaRealNew[kC][kT][kA], * hDEtaMixNew[kC][kT][kA],
             * hDEtaRawFlip[kC][kT][kA], * hDEtaRealFlip[kC][kT][kA];  

        // ----------- loading ntrigg:
        double NTrigg[kC][kT];     
        TH1D * hTriggVSum[kC][kT]; 

        void LoadNtrigg(){
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=0; iptt<fNumPtt; iptt++){
                    NTrigg[ic][iptt]=0;
                    hTriggVSum[ic][iptt]=(TH1D*)fhTriggPt[ic][0][iptt]->Clone(); hTriggVSum[ic][iptt]->Reset();
                    for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++){
                        hTriggVSum[ic][iptt]->Add( fhTriggPt[ic][iv][iptt] );
                        NTrigg[ic][iptt] += fhTriggPt[ic][iv][iptt]->GetEntries();
                    }
                }
            }
            cout << "loading ntrigg done... " << endl;
        }


        // ----------- loading dEta histos:
        void LoadDEta(){

            int rebin_pta=0;
            TString cta;
            MTools * mt = new MTools();

            // ------------------------------
            // loading histos and summing phi/vertex bins
            for(int iptt=1; iptt<fNumPtt; iptt++){
                for(int ipta=0; ipta<fNumPta; ipta++){
                    if( fPTt->At(iptt) < fPTa->At(ipta) )
                        continue;
                    for(int ic=0; ic<fNumCent; ic++){
                        cta = Form("C%.0fT%.0fA%.0f",fCent->At(ic), fPTt->At(iptt), fPTa->At(ipta) );

                        hDEtaRaw[ic][iptt][ipta] = (TH1D*)fhDEtaNearRaw[ic][fVertexSkip][0][iptt][ipta]->Clone("hDEtaRawOLD_"+cta);
                        hDEtaRaw[ic][iptt][ipta]->Reset();
                        hDEtaMix[ic][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][fVertexSkip][0][iptt][ipta]->Clone("hDEtaMixOLD_"+cta);
                        hDEtaMix[ic][iptt][ipta]->Reset();

                        // summing vertex and phiGap
                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++) {
                            for(int ip=0; ip<(fPhiSkip); ip++) {
                                hDEtaRaw[ic][iptt][ipta]->Add( (TH1D*)fhDEtaNearRaw[ic][iv][ip][iptt][ipta] );
                                hDEtaMix[ic][iptt][ipta]->Add( (TH1D*)fhDEtaNearMix[ic][iv][ip][iptt][ipta] );
                            } // Phi
                        } // Vtx
                    } // Cent
                } // PTa
            } // PTt

            // ------------------------------
            // rebinning pTa bins
            // Currently       0.3, 1, 2, 3, 4, 5, 6, 8, 10
            // New should be   2, 3, 4, 6, 8 (specified above)

            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=1; iptt<fNumPtt; iptt++){
                    for(int ipta_new=0; ipta_new<fNumPtaPublished; ipta_new++){
                        if( fPTt->At(iptt) < fPtaPublishedBins[ipta_new] )
                            continue;
                        cta = Form("C%.0fT%.0fA%.0f",fCent->At(ic), fPTt->At(iptt), fPtaPublishedBins[ipta_new] );

                        hDEtaRawNew[ic][iptt][ipta_new] = (TH1D*) hDEtaRaw[ic][iptt][ipta_new]->Clone("hDEtaRaw_"+cta);
                        hDEtaRawNew[ic][iptt][ipta_new]->Reset();

                        hDEtaMixNew[ic][iptt][ipta_new] = (TH1D*) hDEtaMix[ic][iptt][ipta_new]->Clone("hDEtaMix_"+cta);
                        hDEtaMixNew[ic][iptt][ipta_new]->Reset();

                        rebin_pta = 0; // assuming similar bin widths, divide by the number of bins merged
                        for(int ipta=0; ipta<fNumPta; ipta++){
                            if(fPTa->At(ipta+1)<=fPtaPublishedBins[ipta_new] or fPTa->At(ipta)>=fPtaPublishedBins[ipta_new+1])
                                continue;

                            // omit pta 5-6, seems to be problematic
                            if( fPTa->At(ipta)==5 && fPTt->At(iptt)==4) continue;
                            // ---------------------------------
                            cout << fPTt->At(iptt) << "-" << fPTt->At(iptt+1) << endl;
                            cout << "\t" << fPtaPublishedBins[ipta_new] << "-" << fPtaPublishedBins[ipta_new+1] << endl;
                            cout << "\t\t" << fPTa->At(ipta) << "-" <<  fPTa->At(ipta+1) << endl;

                            // ---------------------------------
                            ++rebin_pta;

                            cout << "\t\t\t";
                            hDEtaRaw[ic][iptt][ipta]->Print();
                            cout << "\t\t\t";
                            hDEtaMix[ic][iptt][ipta]->Print();

                            hDEtaRawNew[ic][iptt][ipta_new]->Add( hDEtaRaw[ic][iptt][ipta] );
                            hDEtaMixNew[ic][iptt][ipta_new]->Add( hDEtaMix[ic][iptt][ipta] );

                        } // PTa
                    } // PTa rebinned
                } // PTt
            } // Cent

            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=1; iptt<fNumPtt; iptt++){
                    for(int ipta_new=0; ipta_new<fNumPtaPublished; ipta_new++){
                        if( fPTt->At(iptt) < fPtaPublishedBins[ipta_new] )
                            continue;
                        // norm. to have the value 1 (not the integral)
                        hDEtaMixNew[ic][iptt][ipta_new]->Scale( (2.*1.6) / hDEtaMixNew[ic][iptt][ipta_new]->Integral(), "width" ); 

                        // correct with lowest PTa mixed event, since substantially 
                        // similar to higher PTa bins with much better stat.
                        hDEtaRealNew[ic][iptt][ipta_new]=(TH1D*)hDEtaRawNew[ic][iptt][ipta_new]->Clone();
                        hDEtaRealNew[ic][iptt][ipta_new]->Divide( hDEtaMixNew[ic][iptt][1] );

                        // scale with ntrigg. and correct for rebinning
                        hDEtaRealNew[ic][iptt][ipta_new]->Scale(1./NTrigg[ic][iptt], "width" );
                        hDEtaRawNew[ic][iptt][ipta_new]->Scale(1./NTrigg[ic][iptt], "width" );

//                        hDEtaRawNew[ic][iptt][ipta_new]->Scale(1./float(rebin_pta));
//                        hDEtaMixNew[ic][iptt][ipta_new]->Scale(1./float(rebin_pta));
//                        hDEtaRealNew[ic][iptt][ipta_new]->Scale(1./float(rebin_pta));
                    } // PTa rebinned
                } // PTt
            } // Cent

            cout << "dividing done, now flipping... " << endl;
            // flipping eta:
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=1; iptt<fNumPtt; iptt++){
                    for(int ipta_new=0; ipta_new<fNumPtaPublished; ipta_new++){
                        if( fPTt->At(iptt) < fPtaPublishedBins[ipta_new] )
                            continue;
                        hDEtaRawFlip[ic][iptt][ipta_new] = mt->Flip( hDEtaRawNew[ic][iptt][ipta_new] );
                        hDEtaRealFlip[ic][iptt][ipta_new] = mt->Flip( hDEtaRealNew[ic][iptt][ipta_new] );
                    }
                }
            }
            cout << "rebinning done, starting fitting" << endl;
        }
    
        // ----------- fit rebinned, flipped histos 
        MFit * mfit_gc[5][5][5];    // Gauss+Const.: [Cent][PTt][PTa]
        MFit * mfit_ggc[5][5][5];   // Gen.Gauss+Const.: [Cent][PTt][PTa]
        MFit * mfit_kc[5][5][5];    // Kaplan+Const: [Cent][PTt][PTa]

        void FitRebinned(TString opt="EMRNS") {
            int ifit=0;
            int int_binmin, int_binmax;
            double val, valerr;
            double fitmin=0, fitmax=1.6;
            double int_min=0, int_max=0.4;
            TH1D * htmp;
            
            MTools * mt = new MTools();

            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=1; iptt<fNumPtt; iptt++){
                    hYield[0][ic][iptt] = new TH1D(Form("hYield_GaussConst_C0%dT0%d",ic,iptt), "", fNumPtaPublished, fPtaPublishedBins);
                    hYield[1][ic][iptt] = new TH1D(Form("hYield_GenGaussConst_C0%dT0%d",ic,iptt), "", fNumPtaPublished, fPtaPublishedBins);
                    hYield[2][ic][iptt] = new TH1D(Form("hYield_KaplanConst_C0%dT0%d",ic,iptt), "", fNumPtaPublished, fPtaPublishedBins);

                    hYield_Int[0][ic][iptt] = new TH1D(Form("hYield_fromIntegral_GaussConst_C0%dT0%d",ic,iptt), "", fNumPtaPublished, fPtaPublishedBins);
                    hYield_Int[1][ic][iptt] = new TH1D(Form("hYield_fromIntegral_GenGaussConst_C0%dT0%d",ic,iptt), "", fNumPtaPublished, fPtaPublishedBins);
                    hYield_Int[2][ic][iptt] = new TH1D(Form("hYield_fromIntegral_KaplanConst_C0%dT0%d",ic,iptt), "", fNumPtaPublished, fPtaPublishedBins);

                    for(int ipta=0; ipta<fNumPtaPublished; ipta++){
                        if( fPTt->At(iptt) < fPtaPublishedBins[ipta] )
                            continue;

                        mfit_gc[ic][iptt][ipta]  = new MFit(0,0,hDEtaRealFlip[ic][iptt][ipta], fitmin, fitmax, true);
                        mfit_ggc[ic][iptt][ipta] = new MFit(0,0,hDEtaRealFlip[ic][iptt][ipta], fitmin, fitmax, false);
                        mfit_kc[ic][iptt][ipta]  = new MFit(2,0,hDEtaRealFlip[ic][iptt][ipta], fitmin, fitmax);

                        MFit * mfit_gc_container = new MFit(0,0,hDEtaRealFlip[ic][iptt][ipta], fitmin, fitmax, true);
                        MFit * mfit_kc_container = new MFit(2,0,hDEtaRealFlip[ic][iptt][ipta], fitmin, fitmax);

                        // Gauss + Constant 
                        hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_gc[ic][iptt][ipta]->ffit, opt );
                        hYield[0][ic][iptt]->SetBinContent(ipta+1, mfit_gc[ic][iptt][ipta]->ffit->GetParameter(1) );
                        hYield[0][ic][iptt]->SetBinError(ipta+1, mfit_gc[ic][iptt][ipta]->ffit->GetParError(1) );
                        //mfit_gc_container->ffit->SetParameters( mfit_gc[ic][iptt][ipta]->ffit->GetParameters() );
                        //mfit_gc_container->ffit->SetParameter(0,0);
                        //hYield[0][ic][iptt]->SetBinContent(ipta+1, mfit_gc_container->ffit->Integral(0, 0.6 ) );
                        //hYield[0][ic][iptt]->SetBinError(ipta+1, mfit_gc_container->ffit->IntegralError(0, 0.6) );

                        htmp = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone();
                        int_binmin = htmp->GetXaxis()->FindBin(int_min);
                        int_binmax = htmp->GetXaxis()->FindBin(int_max);
                        mt->subtractConstTH1( htmp,  mfit_gc[ic][iptt][ipta]->ffit->GetParameter(0) );
                        val = htmp->IntegralAndError(int_binmin, int_binmax, valerr);
                        hYield_Int[0][ic][iptt]->SetBinContent(ipta+1, val);
                        hYield_Int[0][ic][iptt]->SetBinError(ipta+1, valerr);

                        // Generalized Gauss + Constant 
                        hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_ggc[ic][iptt][ipta]->ffit, opt );
                        hYield[1][ic][iptt]->SetBinContent(ipta+1,mfit_ggc[ic][iptt][ipta]->ffit->GetParameter(1) );
                        hYield[1][ic][iptt]->SetBinError(ipta+1,mfit_ggc[ic][iptt][ipta]->ffit->GetParError(1) );
                        //mfit_gc_container->ffit->SetParameters( mfit_gc[ic][iptt][ipta]->ffit->GetParameters() );
                        //mfit_gc_container->ffit->SetParameter(0,0);
                        //hYield[1][ic][iptt]->SetBinContent(ipta+1, mfit_gc_container->ffit->Integral(0, 0.6 ) );
                        //hYield[1][ic][iptt]->SetBinError(ipta+1, mfit_gc_container->ffit->IntegralError(0, 0.6) );

                        htmp = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone();
                        mt->subtractConstTH1( htmp,  mfit_ggc[ic][iptt][ipta]->ffit->GetParameter(0) );
                        val = htmp->IntegralAndError(int_binmin, int_binmax, valerr);
                        hYield_Int[1][ic][iptt]->SetBinContent(ipta+1, val);
                        hYield_Int[1][ic][iptt]->SetBinError(ipta+1, valerr);

                        // Kaplan + Constant 
                        hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_kc[ic][iptt][ipta]->ffit, opt );
                        mfit_kc_container->ffit->SetParameters( mfit_kc[ic][iptt][ipta]->ffit->GetParameters() );
                        mfit_kc_container->ffit->SetParameter(0, 0); // set const. to 0 then integrate
                        val = mfit_kc_container->ffit->Integral(0, 0.6);
                        valerr = mfit_kc_container->ffit->IntegralError(0, 0.6);
                        hYield[2][ic][iptt]->SetBinContent(ipta+1, val);
                        hYield[2][ic][iptt]->SetBinError(ipta+1, valerr );

                        htmp = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone();
                        mt->subtractConstTH1( htmp,  mfit_kc[ic][iptt][ipta]->ffit->GetParameter(0) );
                        val = htmp->IntegralAndError(int_binmin, int_binmax, valerr);
                        hYield_Int[2][ic][iptt]->SetBinContent(ipta+1, val);
                        hYield_Int[2][ic][iptt]->SetBinError(ipta+1, valerr);

                        // fix constant of Kaplan from Gaussian fit:
                        // mfit_kc[ic][iptt][ipta]->ffit->FixParameter(0, mfit_gc[ic][iptt][ipta]->ffit->GetParameter(0) );
                        
                    } // PTa
                } // PTt
            } // Cent
        }
        
        // ----------- create output and write histos and fit
        void WriteOutput(TString period, TString AOD, TString RL){
            TString outdir = "../data/processed";

            TString phi = Form("%.1f", fPhi->At(fPhiSkip));
            TString vert= Form("%.1f", fVtx->At(fVertexSkip));
            phi.ReplaceAll(".", "dot"); vert.ReplaceAll(".", "dot");

            if(type == "PbPb") {
                fOutFileName = Form("%s/simplecorr_%s_%sAOD%sRL%s_%s_PhiCut%sVSkip%s.root", outdir.Data(), type.Data(), period.Data(), AOD.Data(), RL.Data(), TrackCut.Data(), phi.Data(), vert.Data());
            }
            if(type == "pp") {
                fOutFileName = Form("%s/simplecorr_%s_%sAOD%s_%s_PhiCut%s.root", outdir.Data(), type.Data(), period.Data(), AOD.Data(), TrackCut.Data(), phi.Data());

            }
            TFile * fOutFile = new TFile( fOutFileName, "RECREATE" );
            fOutFile->cd();
            cout << "Writing output to: " << fOutFileName << endl;

            // write ntrigg and yield
            TString ct, cta;
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=1; iptt<fNumPtt; iptt++){
                    ct = Form("C%.0fT%.0f",fCent->At(ic), fPTt->At(iptt) );

                    hTriggVSum[ic][iptt]->Write("hTrigg_"+ct);

                    hYield[0][ic][iptt]->Write("hYield_GC_"+ct);
                    hYield[1][ic][iptt]->Write("hYield_GGC_"+ct);
                    hYield[2][ic][iptt]->Write("hYield_KC_"+ct);

                    hYield_Int[0][ic][iptt]->Write("hYield_INT_GC_"+ct);
                    hYield_Int[1][ic][iptt]->Write("hYield_INT_GGC_"+ct);
                    hYield_Int[2][ic][iptt]->Write("hYield_INT_KC_"+ct);
                }
            }

            // write eta histograms and fits (only rebinned)
            for(int iptt=1; iptt<fNumPtt; iptt++){
                for(int ipta=0; ipta<fNumPtaPublished; ipta++){
                    if(fPTt->At(iptt) < fPtaPublishedBins[ipta]) 
                        continue;
                    for(int ic=0; ic<fNumCent; ic++){
                        cta = Form("C%.0fT%.0fA%.0f",fCent->At(ic), fPTt->At(iptt), fPtaPublishedBins[ipta]);

                        hDEtaRawNew[ic][iptt][ipta]->Write("hDEtaRaw_"+cta); 
                        hDEtaMixNew[ic][iptt][ipta]->Write("hDEtaMix_"+cta);
                        hDEtaRealNew[ic][iptt][ipta]->Write("hDEtaReal_"+cta);
                        hDEtaRawFlip[ic][iptt][ipta]->Write("hDEtaRawFlip_"+cta);
                        hDEtaRealFlip[ic][iptt][ipta]->Write("hDEtaRealFlip_"+cta);
                        mfit_gc[ic][iptt][ipta]->ffit->Write("ffit_GC_"+cta);
                        mfit_ggc[ic][iptt][ipta]->ffit->Write("ffit_GGC_"+cta);
                        mfit_kc[ic][iptt][ipta]->ffit->Write("ffit_KC_"+cta);
                    }
                }
            }
            fhst->WriteConfigToNew( fOutFile );
            fOutFile->Write();
            fOutFile->Close();
        }
};

void processSimpleCorr() {

    TString jOutDir = "../data/processed/";
    TString jInDir  = "../data/jcorran/";
    
    // info strings related to input file:
    TString InFileName[2]; TString Period[2]; TString AOD[2]; TString RL[2]; TString Comment[2];TString TrackCut[2];
    
    // Loading PbPb
    //
//    InFileName[1] = "JDiHadronCorr_legotrain_allTrigg-CF_PbPb-1815_20151006-1625_runlist_1-LHC10h_AOD86_MgFpMgFm.root"; 
//
    InFileName[1] = "JDiHadronCorr_legotrain_allTrigg-CF_PbPb-1599_20150827-0946_runlist_2-LHC10h_AOD86_MgFpMgFm.root"; //this contains TPCOnly cut
    //InFileName[1] = "JDiHadronCorr_legotrain_allTrigg-CF_PbPb-1815_20151006-1625_runlist_3-LHC10h_AOD86_MgFpMgFm.root"; //this contains TPCOnly cut

    Period[1]="LHC10h"; AOD[1]="86"; RL[1]="2"; Comment[1]=""; TrackCut[1]="TPCOnly";    

    // Loading pp
//    InFileName[0] = "JDiHadronCorr_legotrain_allTrigg-CF_pp-726_20150918-2352-2760GeV_LHC11a_p4_AOD113_noSDD.root";
//    Period[0]="LHC11a"; AOD[0]="113"; RL[0]="0"; Comment[0]="noSDD"; TrackCut[0]="GlobalSDD";    
    InFileName[0] = "JDiHadronCorr_legotrain_allTrigg-CF_pp-727_20150918-2352-2760GeV_LHC11a_p4_AOD113_withSDD.root";
    Period[0]="LHC11a"; AOD[0]="113"; RL[0]="0"; Comment[0]="withSDD"; TrackCut[0]="GlobalSDD";    

    MSimpleCorr * msc[2][5][5];     // [Type][VertexSkip][PhiSkip]

    TString fitopt = "ERNS";
    for(int itype=0; itype<2; itype++){
        msc[itype][0][0] = new MSimpleCorr(jInDir+InFileName[itype], TrackCut[itype], -8, 1.6 );
        msc[itype][0][0]->LoadDEta();
        msc[itype][0][0]->FitRebinned(fitopt);
        msc[itype][0][0]->WriteOutput(Period[itype], AOD[itype], RL[itype]);
        delete msc[itype][0][0];

        msc[itype][0][1] = new MSimpleCorr(jInDir+InFileName[itype], TrackCut[itype], -8, 0.5 );
        msc[itype][0][1]->LoadDEta();
        msc[itype][0][1]->FitRebinned(fitopt);
        msc[itype][0][1]->WriteOutput(Period[itype], AOD[itype], RL[itype]);
        delete msc[itype][0][1];

        msc[itype][0][2] = new MSimpleCorr(jInDir+InFileName[itype], TrackCut[itype], -8, 0.2 );
        msc[itype][0][2]->LoadDEta();
        msc[itype][0][2]->FitRebinned(fitopt);
        msc[itype][0][2]->WriteOutput(Period[itype], AOD[itype], RL[itype]);
        delete msc[itype][0][2];
    }
}
