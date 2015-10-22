#include "../inc/mplot.h"
#include "../inc/mtools.h"
#include "../inc/AliJHistManagerROOT6.cxx"
#include "../inc/loadFilip.h"

// -------------------------------------
// published near-side I_AA, pTt 8-15GeV, C:0-5%
// http://arxiv.org/pdf/1110.0121.pdf
//
TGraphAsymmErrors * publishedIAA(int icase) {
    double xval[] = {3.5, 5.0, 7.0, 9.0};
    double xerrminus[] = {0.5, 1.0, 1.0, 1.0};
    double xerrplus[] = {0.5, 1.0, 1.0, 1.0};
    double yval1[] = {1.4, 1.21, 1.12, 1.25};
    double yerrminus1[] = {0.1118033988749895, 0.09848857801796104, 0.08944271909999159, 0.12041594578792295};
    double yerrplus1[] = {0.1118033988749895, 0.09848857801796104, 0.08944271909999159, 0.12041594578792295};
    int numpoints1 = 4;
    TGraphAsymmErrors * p8088_d1x1y1 = new TGraphAsymmErrors(numpoints1, xval, yval1, xerrminus, xerrplus, yerrminus1, yerrplus1);
    p8088_d1x1y1->SetName("/HepData/8088/d1x1y1");
    p8088_d1x1y1->SetTitle("PRL flat backg.");

    double yval2[] = {1.29, 1.19, 1.12, 1.25};
    double yerrminus2[] ={0.10295630140987, 0.08544003745317531, 0.08944271909999159, 0.12041594578792295};
    double yerrplus2[] = {0.10295630140987, 0.08544003745317531, 0.08944271909999159, 0.12041594578792295};
    int numpoints2 = 4;
    TGraphAsymmErrors * p8088_d1x1y2 = new TGraphAsymmErrors(numpoints2, xval, yval2, xerrminus, xerrplus, yerrminus2, yerrplus2);
    p8088_d1x1y2->SetName("/HepData/8088/d1x1y2");
    p8088_d1x1y2->SetTitle("PRL v2 backg.");

    double yval3[] = {1.25, 1.22, 1.09, 1.3};
    double yerrminus3[] = {0.1140175425099138, 0.09848857801796104, 0.09433981132056604, 0.12727922061357855};
    double yerrplus3[] = {0.1140175425099138, 0.09848857801796104, 0.09433981132056604, 0.12727922061357855};
    int numpoints3 = 4;
    TGraphAsymmErrors * p8088_d1x1y3 = new TGraphAsymmErrors(numpoints3, xval, yval3, xerrminus, xerrplus, yerrminus3, yerrplus3);
    p8088_d1x1y3->SetName("/HepData/8088/d1x1y3");
    p8088_d1x1y3->SetTitle("EtaGap backg.");

    if(icase == 0)
        return p8088_d1x1y1; // 0: flat back.
    if(icase == 1)
        return p8088_d1x1y2; // 1: v2 backg.
    if(icase == 2)
        return p8088_d1x1y3; // 2: EtaGap backg.
}

// -------------------------------------
// gets yield of histogram after subtracting 
// constant from fit
double getYield(TH1 * h1, TF1 * fit){
    MTools * mt = new MTools();
    TH1D * h2 = (TH1D*)h1->Clone();
    mt->subtractConstTH1(h2, fit->GetParameter(0) );
    double yield = h2->Integral(h2->FindBin(1e-6), h2->FindBin(1.1));
    delete h2;
    return yield;
}

void drawSimpleCorr()
{
    // input file (with histogram manager and bins):
    TString inName[2][10]; // [Type][Misc: Phi/Vert]
    TString Types[] = {"pp", "PbPb"};
    cout << "defined types are: " << Types[0] << "\t" << Types[1] << endl;
    TFile * inFile[2][10];

    AliJHistManager * hst[2];
    AliJBin * Cent[2], * Vtx[2], * Eta[2], * Phi[2], * PTt[2], * PTa[2];
    int NumCent[2], NumVtx[2], NumPtt[2], NumPta[2], NumEta[2], NumPhi[2];
    int NSetup[2];
    int NCentrality[] = {1, 5};     // [pp/AA]

    TString jdir = "../data/processed/";
    inName[0][0] = "corr_pp_LHC11aAOD113_GlobalSDD_PhiCut1dot6.root";
    inName[0][1] = "corr_pp_LHC11aAOD113_GlobalSDD_PhiCut0dot5.root";
    inName[0][2] = "corr_pp_LHC11aAOD113_GlobalSDD_PhiCut0dot2.root";
    NSetup[0] = 3;

    inName[1][0] = "corr_PbPb_LHC10hAOD86RL1_TPCOnly_PhiCut1dot6VSkip-8.root";
    inName[1][1] = "corr_PbPb_LHC10hAOD86RL1_TPCOnly_PhiCut0dot5VSkip-8.root";
    inName[1][2] = "corr_PbPb_LHC10hAOD86RL1_TPCOnly_PhiCut0dot2VSkip-8.root";
    NSetup[1] = 3;

    for(int itype=0; itype<2; itype++){
        for(int isetup=0; isetup<NSetup[itype]; isetup++){
            inFile[itype][isetup] = TFile::Open(jdir+inName[itype][isetup]);
        }
    }

    for(int itype=0; itype<2; itype++){
        inFile[itype][0]->cd(); //inFile[itype][0]->ls();
        hst[itype] = new AliJHistManager(Form("hst%d",itype));
        hst[itype]->LoadConfig();
        Cent[itype] = hst[itype]->GetBin("Cent");   NumCent[itype] = Cent[itype]->Size();
        Vtx[itype]  = hst[itype]->GetBin("Vtx");    NumVtx[itype]  = Vtx[itype]->Size();
        PTt[itype]  = hst[itype]->GetBin("PTt");    NumPtt[itype]  = PTt[itype]->Size();
        PTa[itype]  = hst[itype]->GetBin("PTa");    NumPta[itype]  = PTa[itype]->Size();
        Eta[itype]  = hst[itype]->GetBin("EtaGap"); NumEta[itype]  = Eta[itype]->Size();
        Phi[itype]  = hst[itype]->GetBin("PhiGap"); NumPhi[itype]  = Phi[itype]->Size();
    }
    cout << "loading input files done..." << endl;

    // load histos
    const int kF=3; const int kC=5;
    const int kS=4;const int kT=5; const int kA=5;
    TH1D * hDEtaRealFlip[kT][kS][kC][kT][kA];           // [Type][Setup]         [Cent][PTt][PTa]
    TH1D * hDEtaReal[kT][kS][kC][kT][kA];               // [Type][Setup]         [Cent][PTt][PTa]
    TH1D * hDEtaRawFlip[kT][kS][kC][kT][kA];           // [Type][Setup]         [Cent][PTt][PTa]
    TH1D * hDEtaRaw[kT][kS][kC][kT][kA];               // [Type][Setup]         [Cent][PTt][PTa]

//    TH2D * hDPhiReal[

    TH1D * hYield[kT][kS][kF][kC][kT];        // [Type][Setup][FitType][Cent][PTt]
    TH1D * hYield_Int[kT][kS][kF][kC][kT];        // [Type][Setup][FitType][Cent][PTt]
    TH1D * hIAA[kS][kF][kC][kT];                    //       [Setup][FitType][Cent][PTt]
    TH1D * hIAA_Int[kS][kF][kC][kT];                    //       [Setup][FitType][Cent][PTt]
    TF1  * ffit[kT][kS][kF][kC][kT][kA];   // [Type][Setup][FitType][Cent][PTt][Pta]
                
    FilipHistos * hfilip[2];
    hfilip[0] = new FilipHistos("PP");
    hfilip[1] = new FilipHistos("AA");

    double PTtBo[] = {4, 6, 8, 15};
    int NumPttNew = 3;
    double PTaBo[] = {2, 3, 4, 6, 8, 10};
    int NumPtaNew = 5;

    for(int itype=0; itype<2; itype++){
        for(int isetup=0; isetup<NSetup[itype]; isetup++){
            for(int ic=0; ic<NumCent[itype]; ic++){
                for(int iptt=0; iptt<NumPttNew; iptt++){
                // ---------
                    for(int ipta=0; ipta<NumPtaNew; ipta++){
                        if( PTtBo[iptt] < PTaBo[ipta] )
                            continue;

                        hDEtaRealFlip[itype][isetup][ic][iptt][ipta] = (TH1D*)inFile[itype][isetup]->Get(Form("hDEtaRealFlip_C%.0fT%.0fA%.0f", Cent[itype]->At(ic),PTtBo[iptt],PTaBo[ipta]));
                        hDEtaReal[itype][isetup][ic][iptt][ipta] = (TH1D*)inFile[itype][isetup]->Get(Form("hDEtaReal_C%.0fT%.0fA%.0f", Cent[itype]->At(ic),PTtBo[iptt],PTaBo[ipta]));
                        hDEtaRawFlip[itype][isetup][ic][iptt][ipta] = (TH1D*)inFile[itype][isetup]->Get(Form("hDEtaRawFlip_C%.0fT%.0fA%.0f", Cent[itype]->At(ic),PTtBo[iptt],PTaBo[ipta]));
                        hDEtaRaw[itype][isetup][ic][iptt][ipta] = (TH1D*)inFile[itype][isetup]->Get(Form("hDEtaRaw_C%.0fT%.0fA%.0f", Cent[itype]->At(ic),PTtBo[iptt],PTaBo[ipta]));
                        hDEtaRealFlip[itype][isetup][ic][iptt][ipta]->Sumw2(); hDEtaReal[itype][isetup][ic][iptt][ipta]->Sumw2();
                        hDEtaRawFlip[itype][isetup][ic][iptt][ipta]->Sumw2(); hDEtaRaw[itype][isetup][ic][iptt][ipta]->Sumw2();

                        ffit[itype][isetup][0][ic][iptt][ipta] = (TF1*)inFile[itype][isetup]->Get(Form("ffit_GC_C%.0fT%.0fA%.0f",  Cent[itype]->At(ic),PTtBo[iptt],PTaBo[ipta]));
                        ffit[itype][isetup][1][ic][iptt][ipta] = (TF1*)inFile[itype][isetup]->Get(Form("ffit_GGC_C%.0fT%.0fA%.0f", Cent[itype]->At(ic),PTtBo[iptt],PTaBo[ipta]));
                        ffit[itype][isetup][2][ic][iptt][ipta] = (TF1*)inFile[itype][isetup]->Get(Form("ffit_KC_C%.0fT%.0fA%.0f",  Cent[itype]->At(ic),PTtBo[iptt],PTaBo[ipta]));
                    }
                    hYield[itype][isetup][0][ic][iptt] = (TH1D*)inFile[itype][isetup]->Get(Form("hYield_GC_C%.0fT%.0f",  Cent[itype]->At(ic), PTtBo[iptt]));
                    hYield[itype][isetup][1][ic][iptt] = (TH1D*)inFile[itype][isetup]->Get(Form("hYield_GGC_C%.0fT%.0f", Cent[itype]->At(ic), PTtBo[iptt]));
                    hYield[itype][isetup][2][ic][iptt] = (TH1D*)inFile[itype][isetup]->Get(Form("hYield_KC_C%.0fT%.0f",  Cent[itype]->At(ic), PTtBo[iptt]));
                    hYield[itype][isetup][0][ic][iptt]->Sumw2(); 
                    hYield[itype][isetup][1][ic][iptt]->Sumw2(); 
                    hYield[itype][isetup][2][ic][iptt]->Sumw2(); 

                    hYield_Int[itype][isetup][0][ic][iptt] = (TH1D*)inFile[itype][isetup]->Get(Form("hYield_INT_GC_C%.0fT%.0f", Cent[itype]->At(ic), PTtBo[iptt]));
                    hYield_Int[itype][isetup][1][ic][iptt] = (TH1D*)inFile[itype][isetup]->Get(Form("hYield_INT_GGC_C%.0fT%.0f",Cent[itype]->At(ic), PTtBo[iptt]));
                    hYield_Int[itype][isetup][2][ic][iptt] = (TH1D*)inFile[itype][isetup]->Get(Form("hYield_INT_KC_C%.0fT%.0f", Cent[itype]->At(ic), PTtBo[iptt]));
                    hYield_Int[itype][isetup][0][ic][iptt]->Sumw2();
                    hYield_Int[itype][isetup][1][ic][iptt]->Sumw2();
                    hYield_Int[itype][isetup][2][ic][iptt]->Sumw2();
                } // PTt
            } // Cent
        } // Setup
    } // Type

    // do IAA
    for(int ifit=0; ifit<kF; ifit++){ 
        for(int isetup=0; isetup<NSetup[1]; isetup++){
            for(int ic=0; ic<NumCent[1]; ic++){
                for(int iptt=0; iptt<NumPttNew; iptt++){
                    hIAA[isetup][ifit][ic][iptt] = (TH1D*)hYield[1][isetup][ifit][ic][iptt]->Clone();
                    hIAA[isetup][ifit][ic][iptt]->Divide( hYield[0][isetup][ifit][0][iptt] );
                    hIAA[isetup][ifit][ic][iptt]->Print();

                    hIAA_Int[isetup][ifit][ic][iptt] = (TH1D*)hYield_Int[1][isetup][ifit][ic][iptt]->Clone();
                    hIAA_Int[isetup][ifit][ic][iptt]->Divide( hYield_Int[0][isetup][ifit][0][iptt] );
                }
            }
        }
    }
    cout << "loading input done..." << endl;
    
    int iplot = 0; 
    // ___________________________
    //  plot correlation histos
    // ___________________________
    std::vector<TH1*>    hList;  
    std::vector<TString> legList;
    double ymin, ymax, height;
    double const_GC, const_KC, const_GGC;
    for(int itype=0; itype<2; itype++){
        for(int isetup=0; isetup<NSetup[itype]; isetup++){
            for(int ic=0; ic<NumCent[itype]; ic++){
                for(int iptt=0; iptt<NumPttNew; iptt++){
                // ---------
                    for(int ipta=0; ipta<NumPtaNew; ipta++){
                        if( PTtBo[iptt] <= PTaBo[ipta] )
                            continue;

                        MPlot * meta = new MPlot(iplot++, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", false);
                        hDEtaRealFlip[itype][isetup][ic][iptt][ipta]->Rebin(4); hDEtaRealFlip[itype][isetup][ic][iptt][ipta]->Scale(1./4.);
                        hDEtaRawFlip[itype][isetup][ic][iptt][ipta]->Rebin(4); hDEtaRawFlip[itype][isetup][ic][iptt][ipta]->Scale(1./4.);
                        hList   = { hDEtaRealFlip[itype][isetup][ic][iptt][ipta] };
                        legList = {"data"};
                        meta->addHList(hList, legList);
                        // set plot borders from fit instead of histogram
                        const_GC  = ffit[itype][isetup][0][ic][iptt][ipta]->GetParameter(0);
                        const_GGC = ffit[itype][isetup][1][ic][iptt][ipta]->GetParameter(0);
                        const_KC  = ffit[itype][isetup][2][ic][iptt][ipta]->GetParameter(0);
                        TLine * lGC  = new TLine(0, const_GC, 1.6, const_GC );
                        TLine * lGGC = new TLine(0, const_GGC, 1.6, const_GC );
                        TLine * lKC  = new TLine(0, const_KC, 1.6, const_KC );

                        ymax = ffit[itype][isetup][1][ic][iptt][ipta]->Eval(0);
                        ymin = const_KC;
                        height = ymax-ymin;
                        
                        meta->SetLimitsXY(0, 1.6, ymin-height*0.1, ymax+height*0.2);
                       
                        meta->AddInfo( Types[itype] );
                        if(itype==1) meta->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1))); 
                        meta->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
                        meta->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTaBo[ipta], PTaBo[ipta+1]) );
                        meta->Draw();

                        ffit[itype][isetup][0][ic][iptt][ipta]->SetLineColor(kGreen+1);
                        lGC->SetLineColor(kGreen+1);
                        ffit[itype][isetup][1][ic][iptt][ipta]->SetLineColor(kYellow+1);
                        lGGC->SetLineColor(kYellow+1);
                        ffit[itype][isetup][2][ic][iptt][ipta]->SetLineColor(kRed+1);
                        lKC->SetLineColor(kRed+1);

                        meta->AddThisEntry( ffit[itype][isetup][0][ic][iptt][ipta], "Gauss+const. fit", "l");
                        meta->AddThisEntry( (TObject*)0, Form("#chi^{2}/NDF=%.1f/%d",ffit[itype][isetup][0][ic][iptt][ipta]->GetChisquare(),ffit[itype][isetup][0][ic][iptt][ipta]->GetNDF()), "" );
                        meta->AddThisEntry( (TObject*)0, Form("yield=%.2f",ffit[itype][isetup][0][ic][iptt][ipta]->GetParameter(1) ),"" );

                        meta->AddThisEntry( ffit[itype][isetup][1][ic][iptt][ipta], "Gen.Gauss+const. fit", "l");
                        meta->AddThisEntry( (TObject*)0, Form("#chi^{2}/NDF=%.1f/%d",ffit[itype][isetup][1][ic][iptt][ipta]->GetChisquare(),ffit[itype][isetup][1][ic][iptt][ipta]->GetNDF()), "" );
                        meta->AddThisEntry( (TObject*)0, Form("yield=%.2f",ffit[itype][isetup][1][ic][iptt][ipta]->GetParameter(1) ), "" );

                        meta->AddThisEntry( ffit[itype][isetup][2][ic][iptt][ipta], "Kaplan+const. fit", "l");
                        meta->AddThisEntry( (TObject*)0, Form("#chi^{2}/NDF=%.1f/%d",ffit[itype][isetup][2][ic][iptt][ipta]->GetChisquare(),ffit[itype][isetup][2][ic][iptt][ipta]->GetNDF()), "" );

                        meta->DrawThisTF1( (TF1*)ffit[itype][isetup][0][ic][iptt][ipta], "same" );
                        meta->DrawThisTF1( (TF1*)ffit[itype][isetup][1][ic][iptt][ipta], "same" );
                        meta->DrawThisTF1( (TF1*)ffit[itype][isetup][2][ic][iptt][ipta], "same" );
                        meta->DrawThis( lGC, "same" );
                        meta->DrawThis( lKC, "same" );

                        meta->Save(Form("figs/PubYield/PubYied_%s_S0%d_C0%dT0%dA0%d", Types[itype].Data(), isetup, ic, iptt, ipta));


                    } // PTa
                } // PTt
            } // Cent
        } // Setup
    } // Type
                

    // plot IAA (pT, assoc) in published ptt bin
    // ___________________________
    TGraphAsymmErrors * gPubIAA_c[3];
    for(int it=0; it<3; it++){ gPubIAA_c[it] = publishedIAA(it);} // published IAA (all 3 backg. subtraction)

    TString fitnames[] = {"Gauss+C fit", "Gen.Gauss+C fit", "Kaplan+C fit"};
    for(int ifit=0; ifit<kF; ifit++){

        for(int ic=0; ic<1; ic++){
            int iptt= 2;
            // draw yield ------------------
            MPlot * myp1 = new MPlot(iplot++, "p_{T, assoc} [GeV]", "yield", true);
            hList   = { hYield[0][2][ifit][ic][iptt], hYield[1][2][ifit][ic][iptt] };
            legList = { "p-p |#phi|<0.2", "Pb-Pb |#phi|<0.2" };
            myp1->addHList(hList, legList, "l");
            myp1->Draw();
            myp1->AddInfo( fitnames[ifit].Data() );
            myp1->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1))); 
            myp1->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
            myp1->SetRatioLimitsXY(2, 10, 0, 2.5);
            myp1->Save(Form("figs/PubYield/Yield_ppAA_FIT%d_C0%d", ifit, ic));
            // draw IAA --------------------
            MPlot * myp = new MPlot(iplot++, "p_{T, assoc} [GeV]", "I_{AA}", false);
            hList   = { hIAA[1][ifit][ic][iptt], hIAA[2][ifit][ic][iptt], hIAA_Int[1][ifit][ic][iptt], hIAA_Int[2][ifit][ic][iptt] };
            legList = { "#phi<0.5 |v_{z}|<8",    "#phi<0.2 |v_{z}|<8", "INT: #phi<0.5 |v_{z}|<8",    "INT: #phi<0.2 |v_{z}|<8" };
            myp->addHList(hList, legList, "l");
            myp->Draw();
            myp->AddInfo( fitnames[ifit].Data() );
            myp->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1))); 
            myp->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
            myp->AddThisEntry(gPubIAA_c[0], "PRL data (flat)");
            myp->AddThisEntry(gPubIAA_c[1], "PRL data (v_{2})");
            myp->AddThisEntry(gPubIAA_c[2], "PRL data (#eta-gap)");
            myp->DrawThisTGraphAsymmErrors( gPubIAA_c[0], "PE same", 24, 1 );
            myp->DrawThisTGraphAsymmErrors( gPubIAA_c[1], "PE same", 25, 1 );
            myp->DrawThisTGraphAsymmErrors( gPubIAA_c[2], "PE same", 26, 1 );
            myp->SetLimitsXY(3,10,0,2.5);
            myp->Save(Form("figs/PubYield/IAA_FIT%d_C0%d", ifit, ic));
        }
    }

    // plot IAA (delta eta)
    // ___________________________
    // ptt: 0:4-6, 1:6-8, 2:8-15
    // pta: 0:2-3, 1:3-4, 2:4-6, 3:6-8

    // our histos:
    // ptt: 0:4-6, 1:6-8, 2:8-15
    // pta: 0:2-3, 1:3-4, 2:4-6, 3:6-8, 4:8-10

    MTools * mt = new MTools();
    cout << "filip's histos loaded..." << endl << "subtract constant and merge cent" << endl;
    for(int itype=0; itype<2; itype++){
        for(int isetup=0; isetup<NSetup[itype]; isetup++){
            for(int iptt=0; iptt<NumPttNew; iptt++){
            // ---------
                for(int ipta=0; ipta<NumPtaNew-1; ipta++){
                    if( PTtBo[iptt] < PTaBo[ipta] )
                        continue;

                    for(int ic=0; ic<NumCent[itype]; ic++){
                        cout << ic << "\t" << iptt << "\t" << ipta << endl;
                        mt->subtractConstTH1( hDEtaRealFlip[itype][isetup][ic][iptt][ipta], ffit[itype][isetup][0][ic][iptt][ipta]->GetParameter(0) ); // subtract const of Gauss fit
                        mt->subtractConstTH1( hDEtaRawFlip[itype][isetup][ic][iptt][ipta], ffit[itype][isetup][0][ic][iptt][ipta]->GetParameter(0) ); // subtract const of Gauss fit
                    }
                    cout << "constant subtracted" << endl;
                    if(itype==1) {
                        hDEtaRealFlip[itype][isetup][1][iptt][ipta]->Add( hDEtaRealFlip[itype][isetup][0][iptt][ipta] ); // merge 0-5+5-10 = 0-10% for PbPb
                        hDEtaRealFlip[itype][isetup][1][iptt][ipta]->Scale(1./2.);
                        hDEtaRawFlip[itype][isetup][1][iptt][ipta]->Add( hDEtaRawFlip[itype][isetup][0][iptt][ipta] ); // merge 0-5+5-10 = 0-10% for PbPb
                        hDEtaRawFlip[itype][isetup][1][iptt][ipta]->Scale(1./2.);
                    }
                }
            }
        }
    }
    double CentBo[] = {0, 10, 20, 40, 60};
    cout << "Yield (dEta): ours vs. Filip" << endl;

    for(int itype=0; itype<2; itype++){
        for(int isetup=0; isetup<NSetup[itype]; isetup++){
            for(int ic=0; ic<NumCent[itype]; ic++){

                int ic_ours = 0;
                if(itype==1) ic_ours=ic+1;
                if(ic_ours==4) 
                    break;

                for(int iptt=0; iptt<NumPttNew; iptt++){
                // ---------
                    for(int ipta=0; ipta<NumPtaNew; ipta++){
                        if( PTtBo[iptt] <= PTaBo[ipta] )
                            continue;

                        MPlot * metaFilip = new MPlot(iplot++, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", true);
                        hList   = { hfilip[itype]->hDEta[ic][iptt][ipta], hDEtaRawFlip[itype][isetup][ic_ours][iptt][ipta], hDEtaRealFlip[itype][isetup][ic_ours][iptt][ipta] };
                        legList = {"Filip", "Our raw", "Our real"};
                        metaFilip->addHList(hList, legList, "l", "P", "l"); // fasz: why is this not working??? 
                        metaFilip->Draw();
                        if(itype==0) metaFilip->AddInfo( "pp" );
                        if(itype==1) metaFilip->AddInfo( Form("Cent: %.0f-%.0f %%",CentBo[ic], CentBo[ic+1])); 
                        metaFilip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
                        metaFilip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTaBo[ipta], PTaBo[ipta+1]) );
                        metaFilip->SetRatioLimitsXY(0,1.6,-2,3);
                        metaFilip->Save(Form("figs/PubYield/PubYied_Filip_%s_S0%d_C0%dT0%dA0%d", Types[itype].Data(), isetup, ic_ours, iptt, ipta));
                    }
                }
            }
        }
    }
    // draw ratio plot here
    /*
    for(int ic=1; ic<5; ic++){
        for(int ipta=0; ipta<4; ipta++){
            mep[ic][ipta] = new MPlot(iplot++, "|#Delta#eta|", "yield", true, 0.5);
            hList = { hDEtaFlip[0][0][1][ipta], hDEtaFlip[1][ic][1][ipta] };
            legList = {"pp", "PbPb" };
            mep[ic][ipta]->addHList(hList, legList, "PE", "p");
            mep[ic][ipta]->SetLimitsX(0, 0.3);
            mep[ic][ipta]->SetRatioLimitsXY(0,0.3, 0,2);
            mep[ic][ipta]->SetRatioLabel("I_{AA}");

            mep[ic][ipta]->Draw();
            mep[ic][ipta]->Save(Form("figs/PubYield/IAA_DEta_C0%dA0%d", ic, ipta));
        }
    }
    */
    cout << "I_AA (dEta): ours vs. Filip" << endl;
    // draw separate plot here
    // our+Filip:
    // ptt = 0: 4-6, 1:6-8, 2:8-15 Filip
    // pta: 0:2-3, 1:3-4, 2:4-6, 3:6-8
    for(int iptt=0; iptt<NumPttNew; iptt++){
        for(int ipta=0; ipta<NumPtaNew-1; ipta++){ // omit last because of Filip
            if( PTtBo[iptt] <= PTaBo[ipta] ) 
                continue;
            for(int ic=1; ic<5; ic++){
                for(int isetup=0; isetup<NSetup[1]; isetup++){
                    hDEtaRealFlip[1][isetup][ic][iptt][ipta]->Divide( hDEtaRealFlip[0][isetup][0][iptt][ipta] );
                    hDEtaRawFlip[1][isetup][ic][iptt][ipta]->Divide( hDEtaRawFlip[0][isetup][0][iptt][ipta] );
                    //if(ic==1 && iptt==1 && ipta==2){
                    //    hDEtaFlip[1][isetup][ic][iptt][ipta]->Print();
                    //}

                }
                MPlot * miaa = new MPlot(iplot++, "|#Delta#eta|", "I_{AA}", true);
                hList = { hfilip[1]->hIAA[ic-1][iptt][ipta], hDEtaRawFlip[1][1][ic][iptt][ipta], hDEtaRealFlip[1][1][ic][iptt][ipta] };
                legList = {"Filip", "RAW: #phi<0.5 |v_{z}|<8 cm", "REAL: #phi<0.5 |v_{z}|<8 cm"};

                miaa->addHList(hList, legList, "PE", "p");
                miaa->SetLimitsXY(0, 0.28, 0.0, 2.5);
                miaa->SetAppearance( hfilip[1]->hIAA[ic-1][iptt][ipta], 2, 2, 24 );
                miaa->Draw();
                
                miaa->SetRatioLimitsXY(0, 0.28, 0.5, 1.5);
                miaa->AddInfo( Form("Cent: %.0f-%.0f %%",CentBo[ic-1], CentBo[ic])); 
                miaa->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
                miaa->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTaBo[ipta], PTaBo[ipta+1]) );
                miaa->Save(Form("figs/PubYield/IAA_DEta_Filip_C0%dT0%dA0%d", ic, iptt, ipta));
            }
        }
    }

    cout << "I_AA (dEta): reproduced from Filip's histos" << endl;

    TH1D * IAA_Filip[kT][kA];
    for(int iptt=0; iptt<NumPttNew; iptt++){
        for(int ipta=0; ipta<NumPtaNew-1; ipta++){ // omit last because of Filip
            if( PTtBo[iptt] <= PTaBo[ipta] ) 
                continue;
            IAA_Filip[iptt][ipta] = (TH1D*) hfilip[1]->hDEta[0][iptt][ipta]->Clone();
            IAA_Filip[iptt][ipta]->Divide( (TH1D*)hfilip[0]->hDEta[0][iptt][ipta]);
            
            
            MPlot * miaa_filip = new MPlot( iplot++, "|#Delta#eta|", "I_{AA}", true);
            hList = { hfilip[1]->hIAA[0][iptt][ipta], IAA_Filip[iptt][ipta] };
            legList={ "Filip's original IAA", "IAA from Filip's histos" };
            miaa_filip->addHList(hList, legList, "PE", "p");
                miaa_filip->SetLimitsXY(0, 0.28, 0.0, 2.5);
            miaa_filip->Draw();
                miaa_filip->SetRatioLimitsXY(0, 0.28, 0.7, 1.3);
            miaa_filip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
            miaa_filip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTaBo[ipta], PTaBo[ipta+1]) );
            miaa_filip->Save(Form("figs/PubYield/IAA_DEta_FilipVSFilip_C0%dT0%dA0%d", 0, iptt, ipta));
        }
    }


}

