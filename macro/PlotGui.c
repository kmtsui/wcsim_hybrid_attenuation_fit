// example.C

#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>

class PlotMainFrame {
   RQ_OBJECT("PlotMainFrame")
private:
   TGMainFrame         *fMain;
   TRootEmbeddedCanvas *fEcanvas;
   TGPopupMenu *fMenuFile;
   TGMenuBar* fMenuBar;
   TGLayoutHints* fMBItemLayout;
   TGLayoutHints* fMBHelpLayout;
   TGTextEntry* fConfigText;
   TGTextEntry* fOutputText;
   TGNumberEntry* fThreads;
   TGNumberEntry* fSeed;
   std::string fCurDir;
   TFile* fIn;
   std::vector<std::string> histList;
   TGVerticalFrame *fHistframe;
public:
   PlotMainFrame(const TGWindow *p,const char* fName,UInt_t w,UInt_t h);
   virtual ~PlotMainFrame();
   void DoDraw();
   void Config();
   void Output();
   void CloseWindow();
   void RunFit();
   void ListHist();
   void Created() { Emit("Created()"); } //*SIGNAL*
};
class FileHandle {

RQ_OBJECT("FileHandle")

private:
   TGTransientFrame *fMain;
   TGShutter        *fShutter;
   TGLayoutHints    *fLayout;
   const TGPicture  *fDefaultPic;
   TFile* fIn;
   std::vector<std::string> histList;
   std::vector<bool> histAxis;
   std::vector<int> histNbinsX;
   std::vector<double> histXmin;
   std::vector<double> histXmax;
   std::vector<TGNumberEntry*> nBinsEntry;
   std::vector<TGNumberEntry*> xminEntry;
   std::vector<TGNumberEntry*> xmaxEntry;
   std::vector<TGColorSelect*> ColorSel;
   std::vector<bool> histPoly;
   std::vector<std::vector<int>> histPolyOrders;
   std::vector<std::vector<double>> histPolyRanges;
   std::vector<TGTextEntry*> PolyOrderEntry;
   std::vector<TGTextEntry*> PolyRangeEntry;
   TCanvas* cFplot;
   TRootEmbeddedCanvas *fEcanvas;
   TGCheckButton* drawSame;
   TGCheckButton* ignoreFixed;
   TGCheckButton* showPrior;
public:
   FileHandle(const char* fName, const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h);
   ~FileHandle();

   //void AddShutterItem(const char *name, shutterData_t *data);

   // slots
   void CloseWindow();
   void HandleButtons();
   void HandleButtons_xaxis();
   void HandleButtons_poly();
   void SetAxis();
   void ResetAxis();
   void SetPoly();
   void ResetPoly();
   std::vector<double> CalcPol(const double* par, std::vector<double> costh_array, std::vector<int> order, std::vector<double> range);
};
FileHandle::FileHandle(const char* fName, const TGWindow *p, const TGWindow *main,
                         UInt_t w, UInt_t h)
{
   // Create transient frame containing a shutter widget.

   if (strlen(fName)>0)
   {
   	fIn = new TFile(fName);
      if (!fIn->IsOpen()){
         std::cout << "Error, could not open input file: " << fName << std::endl;
         CloseWindow();
         return;
      }
   }

   //cFplot = new TCanvas();

   fMain = new TGTransientFrame(p, main, w, h);
   fMain->Connect("CloseWindow()", "FileHandle", this, "CloseWindow()");
   fMain->DontCallClose(); // to avoid double deletions.

   // use hierarchical cleaning
   fMain->SetCleanup(kDeepCleanup);

   // Create canvas widget
   fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,700,500);
   fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX |
                   kLHintsExpandY, 10,10,10,1));

   TGHorizontalFrame* hframe = new TGHorizontalFrame(fMain,200,40);

   TGVerticalFrame* fHistframe = new TGVerticalFrame(hframe,200,40);
   TList* list = fIn->GetListOfKeys();
   for (int i=0;i<list->GetSize();i++)
   {
      std::string objname = (std::string)list->At(i)->GetName();
      if (objname.find("hist_")==0 && objname.find("_result",objname.size()-7)==objname.size()-7)
      {
         TH1D* hist_result = (TH1D*)fIn->Get(objname.c_str());
         histAxis.push_back(false);
         histNbinsX.push_back(hist_result->GetNbinsX());
         histXmin.push_back(0);
         histXmax.push_back(hist_result->GetNbinsX());
         histPoly.push_back(false);

         objname.erase(0,5); objname.erase(objname.size()-7); 
         std::cout<<"Getting hist: "<<objname<<std::endl;
         histList.push_back(objname);

         TGHorizontalFrame *histframe = new TGHorizontalFrame(fHistframe,200,40);
         TGLabel* lab = new TGLabel(histframe, objname.c_str());
         histframe->AddFrame(lab, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
         TGTextButton *prefit = new TGTextButton(histframe,"Draw", histList.size()-1);
         prefit->Connect("Clicked()","FileHandle",this,"HandleButtons()");
         histframe->AddFrame(prefit, new TGLayoutHints(kLHintsCenterX,
                                                5,5,3,4));
         TGTextButton *postfit = new TGTextButton(histframe,"X-axis", histList.size()-1);
         postfit->Connect("Clicked()","FileHandle",this,"HandleButtons_xaxis()");
         histframe->AddFrame(postfit, new TGLayoutHints(kLHintsCenterX,
                                                5,5,3,4));
         //std::cout<<"ColorSel[histList.size()-1]->GetColor()="<<TColor::GetColor(ColorSel[histList.size()-1]->GetColor())<<std::endl;
         TGTextButton *polysetup = new TGTextButton(histframe,"Poly.", histList.size()-1);
         polysetup->Connect("Clicked()","FileHandle",this,"HandleButtons_poly()");
         histframe->AddFrame(polysetup, new TGLayoutHints(kLHintsCenterX,
                                                5,5,3,4));
         ColorSel.push_back(new TGColorSelect(histframe, TColor::Number2Pixel(histList.size()), 0));
         histframe->AddFrame(ColorSel[histList.size()-1]);
         fHistframe->AddFrame(histframe,new TGLayoutHints(kLHintsCenterX,
                                                5,5,3,4));
      }
   }
   hframe->AddFrame(fHistframe, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));

   nBinsEntry.resize(histList.size());
   xminEntry.resize(histList.size());
   xmaxEntry.resize(histList.size());
   histPolyOrders.resize(histList.size());
   histPolyRanges.resize(histList.size());
   PolyOrderEntry.resize(histList.size());
   PolyRangeEntry.resize(histList.size());

   TGVerticalFrame* fOptframe = new TGVerticalFrame(hframe,200,40);
   TGButtonGroup* optGroup = new TGButtonGroup(fOptframe,"Draw option",kVerticalFrame);
   drawSame = new TGCheckButton(optGroup, new TGHotString("Same"));
   ignoreFixed = new TGCheckButton(optGroup, new TGHotString("Ignore fixed"));
   showPrior = new TGCheckButton(optGroup, new TGHotString("Show Prior"));
   fOptframe->AddFrame(optGroup);
   hframe->AddFrame(fOptframe);

   fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));

   fMain->MapSubwindows();
   fMain->Resize(fMain->GetDefaultSize());

   // position relative to the parent's window
   fMain->CenterOnParent();

   fMain->SetWindowName(fName);

   fMain->MapWindow();
   //gClient->WaitFor(fMain);
}

FileHandle::~FileHandle()
{
   // dtor

   gClient->FreePicture(fDefaultPic);
   fMain->DeleteWindow();  // deletes fMain
}

void FileHandle::CloseWindow()
{
   if (fIn!=nullptr) delete fIn;
   delete this;
}

void FileHandle::HandleButtons()
{
   TGButton *btn = (TGButton *) gTQSender;
   printf("Drawing %s\n", histList[btn->WidgetId()].c_str());

   TCanvas *fCanvas = fEcanvas->GetCanvas();
   fCanvas->cd();

   TH1D* hist_val = (TH1D*)fIn->Get(Form("hist_%s_result",histList[btn->WidgetId()].c_str()));
   TH1D* hist_err = (TH1D*)fIn->Get(Form("hist_%s_error_final",histList[btn->WidgetId()].c_str()));
   TH1D* hist_pre = (TH1D*)fIn->Get(Form("hist_%s_prior",histList[btn->WidgetId()].c_str()));
   TH1D* hist_pre_err = (TH1D*)fIn->Get(Form("hist_%s_error_prior",histList[btn->WidgetId()].c_str()));
   TH1D* hist_result;// = new TH1D("","",hist_val->GetNbinsX(),0,hist_val->GetNbinsX());
   if (!histPoly[btn->WidgetId()])
   {
      if (histAxis[btn->WidgetId()])
         hist_result = new TH1D("","",histNbinsX[btn->WidgetId()],histXmin[btn->WidgetId()],histXmax[btn->WidgetId()]);
      else
      {
         hist_result = (TH1D*)hist_val->Clone();
         hist_result->Reset();
      }
      for (int i=1;i<=hist_val->GetNbinsX();i++)
      {
         if (hist_err->GetBinContent(i)>0||!ignoreFixed->IsDown())
         {
            hist_result->SetBinContent(i,hist_val->GetBinContent(i));
            hist_result->SetBinError(i,hist_err->GetBinContent(i));
         }
      }
   }
   else
   {
      hist_result = new TH1D("","",1000,histPolyRanges[btn->WidgetId()].front(),histPolyRanges[btn->WidgetId()].back());
      std::vector<double> pltpts;
      for (int i=1;i<=hist_result->GetNbinsX();i++)
      {
         pltpts.push_back(hist_result->GetBinCenter(i));
      }
      const int npar = hist_val->GetNbinsX();
      double par_val[npar];
      for (int i=1;i<=npar;i++)
      {
         par_val[i-1] = hist_val->GetBinContent(i);
      }
      std::vector<double> valpts = CalcPol(par_val,pltpts,histPolyOrders[btn->WidgetId()],histPolyRanges[btn->WidgetId()]);
      for (int i=1;i<=hist_result->GetNbinsX();i++)
      {
         hist_result->SetBinContent(i,valpts[i-1]);
      }
   }
   hist_result->SetLineColor(TColor::GetColor(ColorSel[btn->WidgetId()]->GetColor()));
   if (drawSame->IsDown())
      hist_result->Draw("same");
   else
      hist_result->Draw();
   if (showPrior->IsDown())
   {
      for (int i=1;i<=hist_val->GetNbinsX();i++)
      {
         if (hist_pre_err->GetBinContent(i)>0)
         {
            TBox* box = new TBox(hist_result->GetBinLowEdge(i),hist_pre->GetBinContent(i)-hist_pre_err->GetBinContent(i),
                                 hist_result->GetBinLowEdge(i+1),hist_pre->GetBinContent(i)+hist_pre_err->GetBinContent(i));
            box->SetFillColor(856);
            box->Draw("same");
         }
      }
      hist_result->Draw("same");
   }
   fCanvas->Update();
}
void FileHandle::HandleButtons_xaxis()
{
   TGButton *btn = (TGButton *) gTQSender;
   //printf("Shutter button %d\n", btn->WidgetId());

   TGTransientFrame* fXaxis = new TGTransientFrame(gClient->GetRoot(), fMain, 200,40);
   fXaxis->SetCleanup(kDeepCleanup);
   // fXaxis->Connect("CloseWindow()", "FileHandle", this, "DeleteWindow()");
   // fXaxis->DontCallClose(); // to avoid double deletions.
   TGHorizontalFrame* hframe = new TGHorizontalFrame(fXaxis,200,40);
   TGLabel* nBinLabel = new TGLabel(hframe, "No. Bins");
   hframe->AddFrame(nBinLabel, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   nBinsEntry[btn->WidgetId()] = new TGNumberEntry(hframe,histNbinsX[btn->WidgetId()],3,-1,TGNumberFormat::kNESInteger,TGNumberFormat::kNEAPositive);
   hframe->AddFrame(nBinsEntry[btn->WidgetId()], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   TGLabel* xminLabel = new TGLabel(hframe, "Min");
   hframe->AddFrame(xminLabel, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   xminEntry[btn->WidgetId()] = new TGNumberEntry(hframe,histXmin[btn->WidgetId()],3,-1,TGNumberFormat::kNESReal,TGNumberFormat::kNEAAnyNumber);
   hframe->AddFrame(xminEntry[btn->WidgetId()], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   TGLabel* xmaxLabel = new TGLabel(hframe, "Max");
   hframe->AddFrame(xmaxLabel, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   xmaxEntry[btn->WidgetId()] = new TGNumberEntry(hframe,histXmax[btn->WidgetId()],3,-1,TGNumberFormat::kNESReal,TGNumberFormat::kNEAAnyNumber);
   hframe->AddFrame(xmaxEntry[btn->WidgetId()], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   TGTextButton *setaxis = new TGTextButton(hframe,"Set", btn->WidgetId());
   setaxis->Connect("Clicked()","FileHandle",this,"SetAxis()");
   hframe->AddFrame(setaxis, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   TGTextButton *resetaxis = new TGTextButton(hframe,"Reset", btn->WidgetId());
   resetaxis->Connect("Clicked()","FileHandle",this,"ResetAxis()");
   hframe->AddFrame(resetaxis, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   fXaxis->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

   fXaxis->MapSubwindows();
   fXaxis->Resize(fXaxis->GetDefaultSize());

   // position relative to the parent's window
   fXaxis->CenterOnParent();

   fXaxis->SetWindowName(histList[btn->WidgetId()].c_str());

   fXaxis->MapWindow();
}
void FileHandle::HandleButtons_poly()
{
   TGButton *btn = (TGButton *) gTQSender;
   //printf("Shutter button %d\n", btn->WidgetId());

   TGTransientFrame* fPoly = new TGTransientFrame(gClient->GetRoot(), fMain, 200,40);
   fPoly->SetCleanup(kDeepCleanup);

   TGVerticalFrame* vframe = new TGVerticalFrame(fPoly,200,40);
   vframe->AddFrame(new TGLabel(vframe,"Polynomial setup (comma separated)"), new TGLayoutHints(kLHintsLeft,5,5,3,4));

   TGHorizontalFrame* hframe = new TGHorizontalFrame(fPoly,200,40);

   TGLabel* orderLabel = new TGLabel(hframe, "Orders");
   hframe->AddFrame(orderLabel, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   PolyOrderEntry[btn->WidgetId()] = new TGTextEntry(hframe,new TGTextBuffer(5),-1);
   hframe->AddFrame(PolyOrderEntry[btn->WidgetId()], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   TGLabel* rangeLabel = new TGLabel(hframe, "Ranges");
   hframe->AddFrame(rangeLabel, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   PolyRangeEntry[btn->WidgetId()] = new TGTextEntry(hframe,new TGTextBuffer(10),-1);
   hframe->AddFrame(PolyRangeEntry[btn->WidgetId()], new TGLayoutHints(kLHintsCenterX,5,5,3,4));

   TGTextButton *setpoly = new TGTextButton(hframe,"Set", btn->WidgetId());
   setpoly->Connect("Clicked()","FileHandle",this,"SetPoly()");
   hframe->AddFrame(setpoly, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   TGTextButton *resetpoly = new TGTextButton(hframe,"Reset", btn->WidgetId());
   resetpoly->Connect("Clicked()","FileHandle",this,"ResetPoly()");
   hframe->AddFrame(resetpoly, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

   vframe->AddFrame(hframe,new TGLayoutHints(kLHintsCenterX,5,5,3,4));

   fPoly->AddFrame(vframe, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

   fPoly->MapSubwindows();
   fPoly->Resize(fPoly->GetDefaultSize());

   // position relative to the parent's window
   fPoly->CenterOnParent();

   fPoly->SetWindowName(histList[btn->WidgetId()].c_str());

   fPoly->MapWindow();
}
std::vector<double> FileHandle::CalcPol(const double* par, std::vector<double> costh_array, std::vector<int> pol_orders, std::vector<double> pol_range)
{
    // Derive the coefficients for each polynomial
    std::vector<std::vector<double>> pol_coeff;
    std::vector<double> pol_p0, pol_p1;
    int par_index = 0;
    for (int i=0;i<pol_orders.size();i++)
    {
        std::vector<double> coeff;
        if (i==0) // for the first polynomial, the coefficients are unconstrained
        {
            for (int j=0;j<=pol_orders[i];j++)
            {
                coeff.push_back(par[j]);
                par_index++;
            } 
        }
        else // for others, we need to match the 0-th and 1-st order derivatives 
        {
            coeff.push_back(pol_p0[i-1]);
            coeff.push_back(pol_p1[i-1]);
            for (int j=2;j<=pol_orders[i];j++)
            {
                coeff.push_back(par[par_index]);
                par_index++;
            } 
        }
        pol_coeff.emplace_back(coeff);
        // std::cout<<"pol "<<i<<": ";
        // for (int j=0;j<pol_coeff[i].size();j++)
        //     std::cout<<pol_coeff[i][j]<<" ";
        // std::cout<<std::endl;
        double p0 = 0;
        double p1 = 0;
        for (int j=0;j<=pol_orders[i];j++) // store the 0-th and 1-st order derivatives at end-point as boundary conditions
        {
            p0 += coeff[j]*TMath::Power(pol_range[i+1]-pol_range[i],j);
            p1 += coeff[j]*j*TMath::Power(pol_range[i+1]-pol_range[i],j-1);
        } 
        pol_p0.push_back(p0);
        pol_p1.push_back(p1);
    }

    std::vector<double> angular_response;
    for (int k=0;k<costh_array.size();k++) // actually calculate the angular response
    {
        double costh = costh_array[k];
        double val = 0;
        for (int i=0;i<pol_orders.size();i++)
        {
            if (costh>=pol_range[i] && costh<pol_range[i+1])
            {
                //std::cout<<"Using pol"<<i<<std::endl;
                for (int j=0;j<=pol_orders[i];j++)
                {
                    val += pol_coeff[i][j]*TMath::Power(costh-pol_range[i],j);
                }
            }
        }
        //std::cout<<"CalcPol:: costh = "<<costh<<", val = "<<val<<std::endl;
        angular_response.push_back(val);
    }
        

    return angular_response;
}
void FileHandle::SetAxis()
{
   TGButton *btn = (TGButton *) gTQSender;
   printf("Set x-axis for %s\n", histList[btn->WidgetId()].c_str());

   histAxis[btn->WidgetId()] = true;
   histNbinsX[btn->WidgetId()] = (int)nBinsEntry[btn->WidgetId()]->GetIntNumber();
   histXmin[btn->WidgetId()] = xminEntry[btn->WidgetId()]->GetNumber();
   histXmax[btn->WidgetId()] = xmaxEntry[btn->WidgetId()]->GetNumber();
}
void FileHandle::ResetAxis()
{
   TGButton *btn = (TGButton *) gTQSender;
   printf("Reset x-axis for %s\n", histList[btn->WidgetId()].c_str());

   TH1D* hist_val = (TH1D*)fIn->Get(Form("hist_%s_result",histList[btn->WidgetId()].c_str()));
   histAxis[btn->WidgetId()] = false;
   histNbinsX[btn->WidgetId()] = hist_val->GetNbinsX();
   histXmin[btn->WidgetId()] = 0;
   histXmax[btn->WidgetId()] = hist_val->GetNbinsX();
   nBinsEntry[btn->WidgetId()]->SetIntNumber(hist_val->GetNbinsX());
   xminEntry[btn->WidgetId()]->SetNumber(0);
   xmaxEntry[btn->WidgetId()]->SetNumber(hist_val->GetNbinsX());
}
void FileHandle::SetPoly()
{
   TGButton *btn = (TGButton *) gTQSender;
   printf("Set polynomial plot for %s\n", histList[btn->WidgetId()].c_str());

   histPoly[btn->WidgetId()] = true;
   std::string wp = PolyOrderEntry[btn->WidgetId()]->GetText();
   std::stringstream ss1(wp);
   histPolyOrders[btn->WidgetId()].clear();
   for(std::string s; std::getline(ss1, s, ',');)
   {
      histPolyOrders[btn->WidgetId()].push_back(std::stoi(s));
   }
   wp = PolyRangeEntry[btn->WidgetId()]->GetText();
   std::stringstream ss2(wp);
   histPolyRanges[btn->WidgetId()].clear();
   for(std::string s; std::getline(ss2, s, ',');)
   {
      histPolyRanges[btn->WidgetId()].push_back(std::stod(s));
   }
}
void FileHandle::ResetPoly()
{
   TGButton *btn = (TGButton *) gTQSender;
   printf("Reset polynomial plot for %s\n", histList[btn->WidgetId()].c_str());

   histPoly[btn->WidgetId()] = false;
   PolyOrderEntry[btn->WidgetId()]->SetText("");
   histPolyOrders[btn->WidgetId()].clear();
   PolyRangeEntry[btn->WidgetId()]->SetText("");
   histPolyRanges[btn->WidgetId()].clear();
}
PlotMainFrame::PlotMainFrame(const TGWindow *p,const char* fName,UInt_t w,UInt_t h) {
   fCurDir = gSystem->pwd();
   //std::cout<<"fCurDir = "<<fCurDir<<std::endl;
   // Create a main frame
   fMain = new TGMainFrame(p,w,h);

   fMain->Connect("CloseWindow()", "PlotMainFrame", this, "CloseWindow()");

//    // a popup menu
//    fMenuFile = new TGPopupMenu(gClient->GetRoot());
//    // adding menu entries
//    fMenuFile->AddEntry("&Open...",M_FILE_OPEN);
//    fMenuFile->AddEntry("&Close", -1);

//    // menu bar
//    fMenuBar = new TGMenuBar(fMain,100,20,kHorizontalFrame);

//    // menu bar item layout hints
//    fMBItemLayout = new TGLayoutHints(kLHintsTop|kLHintsLeft,0,4,0,0);
//    fMBHelpLayout = new TGLayoutHints(kLHintsTop|kLHintsRight);

//    // adding popup menus
//    fMenuBar->AddPopup("&File", fMenuFile, fMBItemLayout);

   // Create canvas widget
   // fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,200,200);
   // fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX |
   //                 kLHintsExpandY, 10,10,10,1));
   // Create a horizontal frame widget with buttons
   TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);
   //TGVerticalFrame *hframe = new TGVerticalFrame(fMain,200,40);
   // TGTextButton *draw = new TGTextButton(hframe,"&Draw");
   // draw->Connect("Clicked()","PlotMainFrame",this,"DoDraw()");
   // hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX,
   //                                          5,5,3,4));
   TGTextButton *config = new TGTextButton(hframe,"&Open File");
   config->Connect("Clicked()","PlotMainFrame",this,"Config()");
   hframe->AddFrame(config, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
   TGTextButton *exit = new TGTextButton(hframe,"&Exit",
                                "gApplication->Terminate(0)");
   hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
   fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));

   // fConfigText = new TGTextEntry(hframe, new TGTextBuffer(100));
   // fConfigText->SetToolTipText("Input config file");
   // fConfigText->Resize(300, fConfigText->GetDefaultHeight());
   // hframe->AddFrame(fConfigText, new TGLayoutHints(kLHintsTop | kLHintsLeft,
   //                                                     10, 2, 2, 2));

   // TGTextButton *output = new TGTextButton(hframe,"&Output");
   // output->Connect("Clicked()","PlotMainFrame",this,"Output()");
   // hframe->AddFrame(output, new TGLayoutHints(kLHintsCenterX,
   //                                          5,5,3,4));
   // fOutputText = new TGTextEntry(hframe, new TGTextBuffer(100));
   // fOutputText->SetToolTipText("Output file");
   // fOutputText->Resize(300, fOutputText->GetDefaultHeight());
   // hframe->AddFrame(fOutputText, new TGLayoutHints(kLHintsTop | kLHintsLeft,
   //                                                     10, 2, 2, 2));
   // fOutputText->SetText("fitoutput.root");

   // TGLabel* fThreadsLabel = new TGLabel(hframe, "No. Threads");
   // hframe->AddFrame(fThreadsLabel, new TGLayoutHints(kLHintsCenterX,
   //                                          5,5,3,4));
   // fThreads = new TGNumberEntry(hframe,1,3,0,TGNumberFormat::kNESInteger,TGNumberFormat::kNEAPositive);
   // hframe->AddFrame(fThreads, new TGLayoutHints(kLHintsCenterX,
   //                                          5,5,3,4));

   // TGTextButton *runfit = new TGTextButton(hframe,"&Run");
   // runfit->Connect("Clicked()","PlotMainFrame",this,"RunFit()");
   // hframe->AddFrame(runfit, new TGLayoutHints(kLHintsCenterX,
   //                                          5,5,3,4));

   // fHistframe = new TGVerticalFrame(fMain,200,40);
   // fHistframe->SetCleanup(kDeepCleanup);
   // fMain->AddFrame(fHistframe, new TGLayoutHints(kLHintsCenterX,
   //                                           2,2,2,2));

   // Set a name to the main frame
   fMain->SetWindowName("Simple Example");

   // Map all subwindows of main frame
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());

   // Map main frame
   fMain->MapWindow();

   Created();

   if (strlen(fName)>0)
   {
   	fIn = new TFile(fName);
      if (!fIn->IsOpen()){
         std::cout << "Error, could not open input file: " << fName << std::endl;
      }
      else
         new FileHandle(fIn->GetName(), gClient->GetRoot(), fMain, 400, 200);
   }
}
void PlotMainFrame::RunFit() {
   TString appName("optical_fit");
   std::string runCommand(Form("optical_fit -c %s -o %s -n %i \n",fConfigText->GetText(),fOutputText->GetText(),(int)fThreads->GetIntNumber()));
   printf("Run: %s",runCommand.c_str());
   gSystem->ChangeDirectory(fCurDir.c_str());
   gSystem->Exec("pwd");
   gSystem->Exec(runCommand.c_str());
}
void PlotMainFrame::DoDraw() {
   // Draws function graphics in randomly chosen interval
   TF1 *f1 = new TF1("f1","sin(x)/x",0,gRandom->Rndm()*10);
   f1->SetLineWidth(3);
   f1->Draw();
   TCanvas *fCanvas = fEcanvas->GetCanvas();
   fCanvas->cd();
   fCanvas->Update();
}
const char *filetypes[] = { "All files",     "*",
                            "ROOT files",    "*.root",
                            "ROOT macros",   "*.C",
                            "Text files",    "*.[tT][xX][tT]",
                            0,               0 };
void PlotMainFrame::Config(){
    const char *ftypes[] = { "ROOT files",    "*.root",
                             "All files",     "*",
                            0,               0 };
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = ftypes;
    fi.fIniDir    = StrDup(dir);
    printf("fIniDir = %s\n", fi.fIniDir);
    new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
    printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
    dir = fi.fIniDir;
    //fConfigText->SetText(fi.fFilename);

    fIn = new TFile(fi.fFilename);
    if (!fIn->IsOpen()){
         std::cout << "Error, could not open input file: " << fi.fFilename << std::endl;
         return;
    }

   new FileHandle(fIn->GetName(), gClient->GetRoot(), fMain, 400, 200);

   //ListHist();
}
int loop=0;
void PlotMainFrame::ListHist()
{
   if (!fIn->IsOpen()){
      std::cout << "Error, could not open input file" << std::endl;
      return;
   }

   std::cout<<"Getting hist now..."<<std::endl;

   //fHistframe->Clear();
   fMain->RemoveFrame(fHistframe);
   if (fHistframe != nullptr)
      delete fHistframe;
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());

   // Map main frame
   fMain->MapWindow();

   fHistframe = new TGVerticalFrame(fMain,200,40);
   //fHistframe->RemoveAll();
   histList.clear();
   TList* list = fIn->GetListOfKeys();
   for (int i=0;i<list->GetSize();i++)
   {
      std::string objname = (std::string)list->At(i)->GetName();
      if (objname.find("hist_")==0 && objname.find("_result",objname.size()-7)==objname.size()-7)
      {
         objname.erase(0,5); objname.erase(objname.size()-7); 
         std::cout<<"Getting hist: "<<objname<<std::endl;
         objname += Form("%i",loop);
         histList.push_back(objname);

         TGHorizontalFrame *histframe = new TGHorizontalFrame(fHistframe,200,40);
         TGLabel* lab = new TGLabel(histframe, objname.c_str());
         histframe->AddFrame(lab, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
         TGTextButton *prefit = new TGTextButton(histframe,"Prefit");
         prefit->Connect("Clicked()","PlotMainFrame",this,"DoDraw()");
         histframe->AddFrame(prefit, new TGLayoutHints(kLHintsCenterX,
                                                5,5,3,4));
         fHistframe->AddFrame(histframe,new TGLayoutHints(kLHintsCenterX,
                                                5,5,3,4));
      }
   }
   loop++;
   fMain->AddFrame(fHistframe, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));


   // // Set a name to the main frame
   // fMain->SetWindowName("Simple Example");

   // Map all subwindows of main frame
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());

   // Map main frame
   fMain->MapWindow();

   new FileHandle(fIn->GetName(), gClient->GetRoot(), fMain, 400, 200);

}
void PlotMainFrame::Output(){
    const char *ftypes[] = { "ROOT files",    "*.root",
                             "All files",     "*",
                            0,               0 };
    static TString dir(".");
    TGFileInfo fio;
    //fio.fFilename = (char*)"fitoutput.root";
    fio.fFileTypes = ftypes;
    fio.fIniDir    = StrDup(dir);
    printf("fIniDir = %s\n", fio.fIniDir);
    new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fio);
    printf("Save file: %s (dir: %s)\n", fio.fFilename, fio.fIniDir);
    dir = fio.fIniDir;
    fOutputText->SetText(fio.fFilename);
}
void PlotMainFrame::CloseWindow()
{
   // Got close message for this MainFrame. Terminates the application.

   gApplication->Terminate();
}
PlotMainFrame::~PlotMainFrame() {
   // Clean up used widgets: frames, buttons, layout hints
   fMain->Cleanup();
   delete fMain;
}


void PlotGui(const char* fName ="") {
   gStyle->SetOptStat(0);
   // Popup the GUI...
   new PlotMainFrame(gClient->GetRoot(),fName,200,200);
}
