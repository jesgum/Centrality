#pragma once

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TH1.h"
#include "TBox.h"
#include "TPad.h"
#include "TPaletteAxis.h"
#include "TStyle.h"

namespace plot
{

inline void init()
{
  gStyle->SetOptStat(0);
}

// The frame of a pad: the first histogram drawn in it (it defines the axes)
inline TH1* getFrame(TPad* pad)
{
  for (TObject* obj : *pad->GetListOfPrimitives()) {
    if (TH1* h = dynamic_cast<TH1*>(obj))
      return h;
  }
  return nullptr;
}

// All histograms drawn directly in a pad, in drawing order ([0] is the frame)
inline std::vector<TH1*> getHistograms(TPad* pad)
{
  std::vector<TH1*> hists;
  for (TObject* obj : *pad->GetListOfPrimitives()) {
    if (TH1* h = dynamic_cast<TH1*>(obj))
      hists.push_back(h);
  }
  return hists;
}

// Find an object by name (and type) in the pad or any of its sub-pads.
// Returns nullptr if not found; use get() to stop with an error instead.
// (TPad::GetPrimitive does not look inside sub-pads.)
template <typename T = TObject>
inline T* find(TPad* pad, const char* name)
{
  for (TObject* obj : *pad->GetListOfPrimitives()) {
    if (!strcmp(obj->GetName(), name)) {
      if (T* found = dynamic_cast<T*>(obj))
        return found;
    }
    if (obj->InheritsFrom(TPad::Class())) {
      if (T* found = find<T>((TPad*)obj, name))
        return found;
    }
  }
  return nullptr;
}

// Print everything in the canvas (index, class, name), indented per sub-pad.
// The frame of each pad is marked. Boxes also show their x-range, plus the bin
// labels they span if the frame has alphanumeric labels.
inline void printCanvasContents(TPad* c, int depth = 0)
{
  TH1* hFrame = getFrame(c);

  int i = 0;
  for (TObject* obj : *c->GetListOfPrimitives()) {
    printf("%*s%3d  %-10s  %s", 4 * depth, "", i++, obj->ClassName(), obj->GetName());
    if (obj == hFrame)
      printf("  <-- frame of %s", c->GetName());
    if (obj->InheritsFrom(TPad::Class())) {
      printf("\n");
      printCanvasContents((TPad*)obj, depth + 1);
      continue;
    }
    if (obj->IsA() == TBox::Class()) { // exact match: TFrame and TPave also inherit from TBox
      TBox* b = (TBox*)obj;
      printf("  x = [%g, %g]", b->GetX1(), b->GetX2());
      if (hFrame && hFrame->GetXaxis()->GetLabels()) {
        // Box edges sit at bin edges, so step half a bin inwards to find the bins
        int firstBin = hFrame->GetXaxis()->FindBin(b->GetX1() + 0.5 * hFrame->GetXaxis()->GetBinWidth(1));
        int lastBin = hFrame->GetXaxis()->FindBin(b->GetX2() - 0.5 * hFrame->GetXaxis()->GetBinWidth(1));
        printf("  (%s - %s)", hFrame->GetXaxis()->GetBinLabel(firstBin), hFrame->GetXaxis()->GetBinLabel(lastBin));
      }
    }
    printf("\n");
  }
}

// Set the vertical extent of all boxes in the canvas (e.g. after changing the y-range)
inline void setBoxHeights(TCanvas* c, double yMin, double yMax)
{
  for (TObject* obj : *c->GetListOfPrimitives()) {
    if (obj->IsA() == TBox::Class()) { // exact match: TFrame and TPave also inherit from TBox
      TBox* b = (TBox*)obj;
      b->SetY1(yMin);
      b->SetY2(yMax);
    }
  }
  c->Modified();
}

// Remove stats boxes from all histograms in the canvas (and its sub-pads)
inline void removeStats(TPad* pad)
{
  for (TObject* obj : *pad->GetListOfPrimitives()) {
    if (TH1* h = dynamic_cast<TH1*>(obj)) {
      h->SetStats(0);
      if (TObject* st = h->GetListOfFunctions()->FindObject("stats")) {
        h->GetListOfFunctions()->Remove(st);
        delete st;
      }
    } else if (obj->InheritsFrom(TPad::Class())) {
      removeStats((TPad*)obj);
    }
  }
  pad->Modified();
}

// Change the margins of a pad loaded from a file. Setting them directly leaves the
// saved TFrame and colz palette at their old positions (two frames in the output),
// so this also moves the palette next to the new frame and repaints once, which
// updates the saved frame to the new margins.
inline void setMargins(TPad* pad, double left, double right, double bottom, double top)
{
  pad->SetLeftMargin(left);
  pad->SetRightMargin(right);
  pad->SetBottomMargin(bottom);
  pad->SetTopMargin(top);

  for (TH1* h : getHistograms(pad)) {
    if (TPaletteAxis* pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette")) {
      pal->SetX1NDC(1 - right + 0.01);
      pal->SetX2NDC(1 - right + 0.05);
      pal->SetY1NDC(bottom);
      pal->SetY2NDC(1 - top);
    }
  }
  pad->Modified();
  pad->Update();
}

// Set the minimum of all histograms in the pad and its sub-pads
// (for a TH2 drawn with colz this is the lower end of the z-axis)
inline void setMinimum(TPad* pad, double min)
{
  for (TObject* obj : *pad->GetListOfPrimitives()) {
    if (TH1* h = dynamic_cast<TH1*>(obj))
      h->SetMinimum(min);
    else if (obj->InheritsFrom(TPad::Class()))
      setMinimum((TPad*)obj, min);
  }
  pad->Modified();
}

// Set the maximum of all histograms in the pad and its sub-pads
// (for a TH2 drawn with colz this is the upper end of the z-axis)
inline void setMaximum(TPad* pad, double max)
{
  for (TObject* obj : *pad->GetListOfPrimitives()) {
    if (TH1* h = dynamic_cast<TH1*>(obj))
      h->SetMaximum(max);
    else if (obj->InheritsFrom(TPad::Class()))
      setMaximum((TPad*)obj, max);
  }
  pad->Modified();
}

// Remove every TLatex whose text contains `text` from the pad and its sub-pads
// (TLatex objects are all named "TLatex"; the text itself is stored as the title)
inline void removeLatex(TPad* pad, const char* text)
{
  std::vector<TObject*> toRemove;
  for (TObject* obj : *pad->GetListOfPrimitives()) {
    if (obj->InheritsFrom(TPad::Class()))
      removeLatex((TPad*)obj, text);
    else if (obj->InheritsFrom(TLatex::Class()) && TString(obj->GetTitle()).Contains(text))
      toRemove.push_back(obj);
  }
  for (TObject* obj : toRemove) {
    pad->GetListOfPrimitives()->Remove(obj);
    delete obj;
  }
  pad->Modified();
}

// Text sizes that are the same in every figure, whatever the canvas or pad size.
// A precision-2 font (e.g. 42) takes SetTextSize as a fraction of the smaller side
// of the pad it is drawn in, so the same value gives different sizes in different
// canvases. ROOT also scales every canvas to a ~20 cm wide PDF page, so a fixed size
// in pixels still differs between canvases of different width. What stays the same
// in the PDFs (and in the note, with all figures at the same width) is the text
// size as a fraction of the canvas width.

// Convert a precision-2 text size in `pad` to a fraction of the canvas width
inline double toCanvasWidthFraction(TPad* pad, double size)
{
  const double padWidth = pad->GetWw() * pad->GetAbsWNDC();
  const double padHeight = pad->GetWh() * pad->GetAbsHNDC();
  return size * std::min(padWidth, padHeight) / pad->GetCanvas()->GetWw();
}

// Set the text size to `fraction` of the canvas width, using `font` (e.g. 4 for
// Helvetica as in 42, 6 for bold as in 62) with pixel precision (43, 63)
inline void setTextSizeOfCanvasWidth(TAttText* text, TCanvas* c, double fraction, int font = 4)
{
  text->SetTextFont(10 * font + 3);
  text->SetTextSize(fraction * c->GetWw());
}

// Get any object by name (and type) from the canvas, searching all sub-pads.
// Stops the macro with a clear message, plus the canvas contents, if it is
// not there, instead of returning nullptr and crashing later.
//   TH1* h      = plot::get<TH1>(c, "hCoT0_Zoom");
//   TLegend* l  = plot::get<TLegend>(c, "TPave");
//   TPad* pad   = plot::get<TPad>(c, "padUp2");
template <typename T = TObject>
inline T* get(TPad* pad, const char* name)
{
  if (T* found = find<T>(pad, name))
    return found;
  printf("plot::get: no %s named \"%s\" in %s. Contents:\n", T::Class_Name(), name, pad->GetName());
  printCanvasContents(pad);
  throw std::runtime_error(Form("plot::get: \"%s\" not found", name));
}

// Get the frame (first histogram drawn) of the named pad, e.g. plot::getFrame(c, "padUp2")
inline TH1* getFrame(TPad* c, const char* padName)
{
  TH1* frame = getFrame(get<TPad>(c, padName));
  if (!frame)
    throw std::runtime_error(Form("plot::getFrame: pad \"%s\" has no histogram", padName));
  return frame;
}

// Draw a horizontal line at y (default 1) across the visible x-range of the pad's
// frame, e.g. at unity in a ratio pad. Call it after changing the x-range. By default
// the line goes behind everything except the frame histogram itself; pass
// behind = false to draw it on top instead.
inline TLine* drawHorizontalLine(TPad* pad, double y = 1., bool behind = true)
{
  TH1* frame = getFrame(pad);
  if (!frame) {
    throw std::runtime_error(Form("plot::drawHorizontalLine: pad \"%s\" has no histogram", pad->GetName()));
  }
  const TAxis* axis = frame->GetXaxis();
  TLine* line = new TLine(axis->GetBinLowEdge(axis->GetFirst()), y, axis->GetBinUpEdge(axis->GetLast()), y);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kGray + 2);
  if (behind) {
    pad->GetListOfPrimitives()->AddAfter(frame, line);
  } else {
    pad->cd();
    line->Draw();
  }
  pad->Modified();
  return line;
}

} // namespace plot
