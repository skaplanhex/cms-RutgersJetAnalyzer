#include "TStyle.h"

void setEXOStyle() {

  TStyle *exoStyle = new TStyle("exoStyle","Style for EXOTICA");

  // For the canvas:
  exoStyle->SetCanvasBorderMode(0);
  exoStyle->SetCanvasColor(kWhite);
  exoStyle->SetCanvasDefH(600); //Height of canvas
  exoStyle->SetCanvasDefW(600); //Width of canvas
  exoStyle->SetCanvasDefX(0);   //POsition on screen
  exoStyle->SetCanvasDefY(0);

  // For the Pad:
  exoStyle->SetPadBorderMode(0);
  // exoStyle->SetPadBorderSize(Width_t size = 1);
  exoStyle->SetPadColor(kWhite);
  exoStyle->SetPadGridX(0);
  exoStyle->SetPadGridY(0);
  exoStyle->SetGridColor(0);
  exoStyle->SetGridStyle(3);
  exoStyle->SetGridWidth(1);

  // For the frame:
  exoStyle->SetFrameBorderMode(0);
  exoStyle->SetFrameBorderSize(1);
  exoStyle->SetFrameFillColor(0);
  exoStyle->SetFrameFillStyle(0);
  exoStyle->SetFrameLineColor(1);
  exoStyle->SetFrameLineStyle(1);
  exoStyle->SetFrameLineWidth(1);

  // For the histo:
  // exoStyle->SetHistFillColor(1);
  // exoStyle->SetHistFillStyle(0);
  exoStyle->SetHistLineColor(1);
  exoStyle->SetHistLineStyle(0);
  exoStyle->SetHistLineWidth(1);
  // exoStyle->SetLegoInnerR(Float_t rad = 0.5);
  // exoStyle->SetNumberContours(Int_t number = 20);

  exoStyle->SetEndErrorSize(2);
  // exoStyle->SetErrorMarker(20);
  // exoStyle->SetErrorX(0.);

  exoStyle->SetMarkerStyle(20);

  //For the fit/function:
  exoStyle->SetOptFit(1);
  exoStyle->SetFitFormat("5.4g");
  exoStyle->SetFuncColor(2);
  exoStyle->SetFuncStyle(1);
  exoStyle->SetFuncWidth(1);

  // For the date:
  exoStyle->SetOptDate(0);
  // exoStyle->SetDateX(Float_t x = 0.01);
  // exoStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  exoStyle->SetOptFile(0);
  exoStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  exoStyle->SetStatColor(kWhite);
  exoStyle->SetStatFont(42);
  exoStyle->SetStatFontSize(0.025);
  exoStyle->SetStatTextColor(1);
  exoStyle->SetStatFormat("6.4g");
  exoStyle->SetStatBorderSize(1);
  exoStyle->SetStatH(0.1);
  exoStyle->SetStatW(0.15);
  // exoStyle->SetStatStyle(Style_t style = 1001);
  // exoStyle->SetStatX(Float_t x = 0);
  // exoStyle->SetStatY(Float_t y = 0);

  // Margins:
  exoStyle->SetPadTopMargin(0.06);
  exoStyle->SetPadBottomMargin(0.14);
  exoStyle->SetPadLeftMargin(0.15);
  exoStyle->SetPadRightMargin(0.05);

  // For the Global title:
  exoStyle->SetOptTitle(0);
  exoStyle->SetTitleFont(42);
  exoStyle->SetTitleColor(1);
  exoStyle->SetTitleTextColor(1);
  exoStyle->SetTitleFillColor(10);
  exoStyle->SetTitleFontSize(0.05);
  // exoStyle->SetTitleH(0); // Set the height of the title box
  // exoStyle->SetTitleW(0); // Set the width of the title box
  // exoStyle->SetTitleX(0); // Set the position of the title box
  // exoStyle->SetTitleY(0.985); // Set the position of the title box
  // exoStyle->SetTitleStyle(Style_t style = 1001);
  // exoStyle->SetTitleBorderSize(2);

  // For the axis titles:
  exoStyle->SetTitleColor(1, "XYZ");
  exoStyle->SetTitleFont(42, "XYZ");
  exoStyle->SetTitleSize(0.06, "XYZ");
  // exoStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // exoStyle->SetTitleYSize(Float_t size = 0.02);
  exoStyle->SetTitleXOffset(0.9);
  exoStyle->SetTitleYOffset(1.25);
  // exoStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  exoStyle->SetLabelColor(1, "XYZ");
  exoStyle->SetLabelFont(42, "XYZ");
  exoStyle->SetLabelOffset(0.007, "XYZ");
  exoStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  exoStyle->SetAxisColor(1, "XYZ");
  exoStyle->SetStripDecimals(kTRUE);
  exoStyle->SetTickLength(0.03, "XYZ");
  exoStyle->SetNdivisions(510, "XYZ");
  exoStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  exoStyle->SetPadTickY(1);

  // Change for log plots:
  exoStyle->SetOptLogx(0);
  exoStyle->SetOptLogy(0);
  exoStyle->SetOptLogz(0);

  // Postscript options:
  exoStyle->SetPaperSize(20.,20.);
  // exoStyle->SetLineScalePS(Float_t scale = 3);
  // exoStyle->SetLineStyleString(Int_t i, const char* text);
  // exoStyle->SetHeaderPS(const char* header);
  // exoStyle->SetTitlePS(const char* pstitle);

  // exoStyle->SetBarOffset(Float_t baroff = 0.5);
  // exoStyle->SetBarWidth(Float_t barwidth = 0.5);
  // exoStyle->SetPaintTextFormat(const char* format = "g");
  // exoStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // exoStyle->SetTimeOffset(Double_t toffset);
  // exoStyle->SetHistMinimumZero(kTRUE);

  exoStyle->cd();

}
