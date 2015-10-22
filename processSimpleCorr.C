/* 
 * This macro compares our results to published
 * http://hepdata.cedar.ac.uk/view/ins930312
 * without correction in vertex bins, without mixed event correction
 * with TPCOnly cut, in pTt 8-15, pTa 3-4, 4-6, 6-8, 8-10 GeV
 */

#include "../inc/AliJHistManagerROOT6.cxx"
#include "../inc/mfit.h"
#include "../inc/mtools.h"

// -------------------------------------
// Simplified correlation analysis class
//
class MSimpleCorr
{
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

            TFile * InFile = TFile::Open(inName);
            InFile->cd(); //InFile->ls();

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

            cout << "initialization done..." << endl;
            LoadNtrigg();
        }

        // output containers
        TH1D * hDEtaRaw[5][5][5];       // [Cent][PTt][Pta]
        TH1D * hDEtaMix[5][5][5];       // [Cent][PTt][Pta]
        TH1D * hDEtaReal[5][5][5];      // [Cent][PTt][Pta]

        // output histos with new pta binning
        TH1D * hYield[3][5][5];         // [FitType][Cent][PTt]
        TH1D * hYield_Int[3][5][5];     // [FitType][Cent][PTt] (yield calculated by integrating histograms. Constans estimated from the different fits.)
        TH1D * hDEtaRawNew[5][5][5];    // [Cent][PTt][Pta] (rebinned for published PTa bins
        TH1D * hDEtaRealNew[5][5][5];   // [Cent][PTt][Pta] (rebinned for published PTa bins
        TH1D * hDEtaMixNew[5][5][5];    // [Cent][PTt][Pta] (rebinned for published PTa bins
        TH1D * hDEtaRawFlip[5][5][5];   // [Cent][PTt][Pta] (rebinned for published PTa bins
        TH1D * hDEtaRealFlip[5][5][5];  // [Cent][PTt][Pta] (rebinned for published PTa bins

        // ----------- loading ntrigg:
        double NTrigg[5][5];     // [Cent][PTt]
        TH1D * hTriggVSum[5][5]; // [Cent][PTt]

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

        // ----------- print bin content to check for zeros:
        void CheckBins() {
            int nbins_out;
            double x,y,yerr;
            int iptt = fPTt->GetBin(8.); 
            int ic = 0;

            for(int ipta=0; ipta<fNumPta; ipta++){
                ofstream out;
                out.open(Form("figs/PubYield/%s_%d_%d.txt", type.Data(),ic,ipta));
                for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++) {
                    for(int ip=0; ip<(fPhiSkip); ip++) {
                        // QA of bin content:
                        out << "##" << endl;
                        out << "iv = " << iv << "  ip = " <<  ip << endl;
                        nbins_out = fhDEtaNearRaw[ic][iv][ip][iptt][ipta]->GetNbinsX();
                        for(int ib=0; ib<nbins_out; ib++){
                            x = fhDEtaNearRaw[ic][iv][ip][iptt][ipta]->GetBinCenter(ib+1);
                            y = fhDEtaNearRaw[ic][iv][ip][iptt][ipta]->GetBinContent(ib+1);
                            yerr = fhDEtaNearRaw[ic][iv][ip][iptt][ipta]->GetBinError(ib+1);
                            out << "\t\t" << ib+1 << "\t" << x << "\t" << y << "\t" << yerr/sqrt(y) << endl;
                        }
                    }
                }
                out.close();
            }
        }

        // ----------- loading dEta histos:
        void LoadDEta(){

            MTools * mt = new MTools();
            for(int iptt=1; iptt<fNumPtt; iptt++){
                for(int ipta=0; ipta<fNumPta; ipta++){
                    if( fPTt->At(iptt) < fPTa->At(ipta) )
                        continue;
                    for(int ic=0; ic<fNumCent; ic++){
                        hDEtaRaw[ic][iptt][ipta] = (TH1D*)fhDEtaNearRaw[ic][fVertexSkip][0][0][0]->Clone();
                        hDEtaRaw[ic][iptt][ipta]->Reset();
                        hDEtaMix[ic][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][fVertexSkip][0][0][0]->Clone();
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

            cout << "summing vertex+phi done... " << endl;
            // rebinning pTa bins
            // Currently       0.3, 1, 2, 3, 4, 5, 6, 8, 10
            // New should be   2, 3, 4, 6, 8 (specified above)
            
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=1; iptt<fNumPtt; iptt++){
                    for(int ipta_new=0; ipta_new<fNumPtaPublished; ipta_new++){
                        if( fPTt->At(iptt) < fPtaPublishedBins[ipta_new] )
                            continue;

                        hDEtaRawNew[ic][iptt][ipta_new] = (TH1D*) hDEtaRaw[ic][iptt][ipta_new]->Clone();
                        hDEtaRawNew[ic][iptt][ipta_new]->Reset();

                        hDEtaMixNew[ic][iptt][ipta_new] = (TH1D*) hDEtaMix[ic][iptt][ipta_new]->Clone();
                        hDEtaMixNew[ic][iptt][ipta_new]->Reset();

                        int rebin_pta = 0; // assuming similar bin widths, divide by the number of bins merged
                        for(int ipta=0; ipta<fNumPta; ipta++){
                            if(fPTa->At(ipta+1)<=fPtaPublishedBins[ipta_new] or fPTa->At(ipta)>=fPtaPublishedBins[ipta_new+1])
                                continue;

                            // omit pta 5-6, seems to be problematic
                            if( fPTa->At(ipta)==5) continue;
                            // ---------------------------------
                            cout << fPTt->At(iptt) << "-" << fPTt->At(iptt+1) << endl;
                            cout << "\t" << fPtaPublishedBins[ipta_new] << "-" << fPtaPublishedBins[ipta_new+1] << endl;
                            cout << "\t\t" << fPTa->At(ipta) << "-" <<  fPTa->At(ipta+1) << endl;
                            cout << "\t\t\t";
                            hDEtaRaw[ic][iptt][ipta]->Print();
                            // ---------------------------------
                            rebin_pta++;
                            if( rebin_pta > 1) cout << "bin should be normalized" << endl;
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
                        ipta_new==0 ? hDEtaRealNew[ic][iptt][ipta_new]->Divide( hDEtaMixNew[ic][iptt][0] ) : hDEtaRealNew[ic][iptt][ipta_new]->Divide( hDEtaMixNew[ic][iptt][1] ); // dividing with pta: 3-4 (except in bin 2-3

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
        TH1D * htmp;
        int ifit=0;

        void FitRebinned() {
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

                        mfit_gc[ic][iptt][ipta]  = new MFit(ifit++,0,0,hDEtaRealFlip[ic][iptt][ipta], 0, 1.6, true);
                        mfit_ggc[ic][iptt][ipta] = new MFit(ifit++,0,0,hDEtaRealFlip[ic][iptt][ipta], 0, 1.6, false);
                        mfit_kc[ic][iptt][ipta]  = new MFit(ifit++,3,0,hDEtaRealFlip[ic][iptt][ipta], 0, 1.6);

                        MFit * mfit_gc_container = new MFit(ifit++,0,0,hDEtaRealFlip[ic][iptt][ipta], 0, 1.6, true);
                        MFit * mfit_kc_container = new MFit(ifit++,3,0,hDEtaRealFlip[ic][iptt][ipta], 0, 1.6);

                        // Gauss + Constant 
                        hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_gc[ic][iptt][ipta]->ffit, "EMRNS" );
                        hYield[0][ic][iptt]->SetBinContent(ipta+1, mfit_gc[ic][iptt][ipta]->ffit->GetParameter(1) );
                        hYield[0][ic][iptt]->SetBinError(ipta+1, mfit_gc[ic][iptt][ipta]->ffit->GetParError(1) );

                        htmp = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone();
                        mt->subtractConstTH1( htmp,  mfit_gc[ic][iptt][ipta]->ffit->GetParameter(0) );
                        hYield_Int[0][ic][iptt]->SetBinContent(ipta+1, htmp->Integral());

                        // Generalized Gauss + Constant 
                        hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_ggc[ic][iptt][ipta]->ffit, "EMRNS" );
                        hYield[1][ic][iptt]->SetBinContent(ipta+1, mfit_ggc[ic][iptt][ipta]->ffit->GetParameter(1) );
                        hYield[1][ic][iptt]->SetBinError(ipta+1, mfit_ggc[ic][iptt][ipta]->ffit->GetParError(1) );

                        htmp = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone();
                        mt->subtractConstTH1( htmp,  mfit_ggc[ic][iptt][ipta]->ffit->GetParameter(0) );
                        hYield_Int[1][ic][iptt]->SetBinContent(ipta+1, htmp->Integral());

                        // Kaplan + Constant 
                        hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_kc[ic][iptt][ipta]->ffit, "EMRNS" );
                        mfit_kc_container->ffit->SetParameters( mfit_kc[ic][iptt][ipta]->ffit->GetParameters() );
                        mfit_kc_container->ffit->SetParameter(0, 0); // set const. to 0 then integrate
                        double val = mfit_kc_container->ffit->Integral(1e-2, 1.1);
                        double valerr = mfit_kc_container->ffit->IntegralError(1e-2, 1.1);
                        hYield[2][ic][iptt]->SetBinContent(ipta+1, val);
                        hYield[2][ic][iptt]->SetBinError(ipta+1, valerr );

                        htmp = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone();
                        mt->subtractConstTH1( htmp,  mfit_kc[ic][iptt][ipta]->ffit->GetParameter(0) );
                        hYield_Int[2][ic][iptt]->SetBinContent(ipta+1, htmp->Integral());

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
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=1; iptt<fNumPtt; iptt++){
                    hTriggVSum[ic][iptt]->Write(Form("hTrigg_C0%dT0%.0f",   ic, fPTt->At(iptt)));

                    hYield[0][ic][iptt]->Write(Form("hYield_GC_C0%dT0%.0f", ic, fPTt->At(iptt)));
                    hYield[1][ic][iptt]->Write(Form("hYield_GGC_C0%dT0%.0f", ic, fPTt->At(iptt)));
                    hYield[2][ic][iptt]->Write(Form("hYield_KC_C0%dT0%.0f", ic, fPTt->At(iptt)));

                    hYield_Int[0][ic][iptt]->Write(Form("hYield_INT_GC_C0%dT0%.0f", ic, fPTt->At(iptt)));
                    hYield_Int[1][ic][iptt]->Write(Form("hYield_INT_GGC_C0%dT0%.0f", ic, fPTt->At(iptt)));
                    hYield_Int[2][ic][iptt]->Write(Form("hYield_INT_KC_C0%dT0%.0f", ic, fPTt->At(iptt)));
                }
            }

            // write eta histograms and fits (only rebinned)
            TString cta = "";
            for(int iptt=1; iptt<fNumPtt; iptt++){
                for(int ipta=0; ipta<fNumPtaPublished; ipta++){
                    if(fPTt->At(iptt) < fPtaPublishedBins[ipta]) 
                        continue;
                    for(int ic=0; ic<fNumCent; ic++){
                        cta = Form("C0%dT0%.fA0%.0f",ic, fPTt->At(iptt), fPtaPublishedBins[ipta]);

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
    InFileName[1] = "JDiHadronCorr_legotrain_allTrigg-CF_PbPb-1815_20151006-1625_runlist_1-LHC10h_AOD86_MgFpMgFm.root"; //this contains TPCOnly cut
    //InFileName[1] = "JDiHadronCorr_legotrain_allTrigg-CF_PbPb-1599_20150827-0946_runlist_2-LHC10h_AOD86_MgFpMgFm.root"; //this contains TPCOnly cut
    //InFileName[1] = "JDiHadronCorr_legotrain_allTrigg-CF_PbPb-1815_20151006-1625_runlist_3-LHC10h_AOD86_MgFpMgFm.root"; //this contains TPCOnly cut

    Period[1]="LHC10h"; AOD[1]="86"; RL[1]="1"; Comment[1]=""; TrackCut[1]="TPCOnly";    

    // Loading pp
    //InFileName[0] = "JDiHadronCorr_legotrain_allTrigg-CF_pp-726_20150918-2352-2760GeV_LHC11a_p4_AOD113_noSDD.root";
    //Period[0]="LHC11a"; AOD[0]="113"; RL[0]="0"; Comment[0]="noSDD"; TrackCut[0]="GlobalSDD";    
    InFileName[0] = "JDiHadronCorr_legotrain_allTrigg-CF_pp-727_20150918-2352-2760GeV_LHC11a_p4_AOD113_withSDD.root";
    Period[0]="LHC11a"; AOD[0]="113"; RL[0]="0"; Comment[0]="withSDD"; TrackCut[0]="GlobalSDD";    

    MSimpleCorr * msc[2][5][5];     // [Type][VertexSkip][PhiSkip]

    for(int itype=0; itype<2; itype++){
        msc[itype][0][0] = new MSimpleCorr(jInDir+InFileName[itype], TrackCut[itype], -8, 1.6 );
        msc[itype][0][0]->LoadDEta();
        msc[itype][0][0]->FitRebinned();
        msc[itype][0][0]->WriteOutput(Period[itype], AOD[itype], RL[itype]);

        msc[itype][0][1] = new MSimpleCorr(jInDir+InFileName[itype], TrackCut[itype], -8, 0.5 );
        msc[itype][0][1]->LoadDEta();
        msc[itype][0][1]->FitRebinned();
        msc[itype][0][1]->WriteOutput(Period[itype], AOD[itype], RL[itype]);

        msc[itype][0][2] = new MSimpleCorr(jInDir+InFileName[itype], TrackCut[itype], -8, 0.2 );
        msc[itype][0][2]->LoadDEta();
        msc[itype][0][2]->FitRebinned();
        msc[itype][0][2]->WriteOutput(Period[itype], AOD[itype], RL[itype]);

        delete msc[itype][0][0];
        delete msc[itype][0][1];
        delete msc[itype][0][2];
    }
}
