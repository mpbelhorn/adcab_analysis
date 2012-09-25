#include "Fit3D.h"

Fit3D::Fit3D(
    TString input_ntuple_file,
    TString analysis_name,
    TString x_axis_label,
    TString x_axis_unit,
    double min_x_bin_edge,
    double max_x_bin_edge,
    TString y_axis_label,
    TString y_axis_unit,
    double min_y_bin_edge,
    double max_y_bin_edge,
    TString z_axis_label,
    TString z_axis_unit,
    double min_z_bin_edge,
    double max_z_bin_edge)
  : DileptonEvents(input_ntuple_file, analysis_name)
{
  std::cout << "Initiliazing Fit3D class." << std::endl;
  
  // Remove information about moving files. Really, all RooMsgService streams
  //   should be sent to a debug log.
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);

  // Define the data set elements and add them to the dataset.
  x_variable_ = new RooRealVar(
      "x_variable", x_axis_label, min_x_bin_edge, max_x_bin_edge, x_axis_unit);
  y_variable_ = new RooRealVar(
      "y_variable", y_axis_label, min_y_bin_edge, max_y_bin_edge, y_axis_unit);
  z_variable_ = new RooRealVar(
      "z_variable", z_axis_label, min_z_bin_edge, max_z_bin_edge, z_axis_unit);
  x_variable_->setBins(30);
  y_variable_->setBins(30);
  z_variable_->setBins(30);
  std::cout << "Adding new columns to dataset." << std::endl;
  RooArgSet analysis_variables(*x_variable_, *y_variable_, *z_variable_);
  data_set_->printArgs(cout);
  std::cout << std::endl;
  recreateDataSet(analysis_variables);
  data_set_->printArgs(cout);
  std::cout << std::endl;
  
  // Construct the histograms.
  TString name;
  name.Clear();
  histograms_.clear();

  for (int i = 0; i < tags_.species_.size(); ++i) {
    HistogramListOverEventSigns event_sign_histograms;
    event_sign_histograms.clear();
    for (int j = 0; j < tags_.signs_.size(); ++j) {
      HistogramListOverComponents component_histograms;
      component_histograms.clear();
      for (int k = 0; k < tags_.components_.size(); ++k) {
        name.Append(tags_.species_[i]);
        name.Append("_");
        name.Append(tags_.signs_[j]);
        name.Append("_");
        name.Append(tags_.components_[k]);
        TH3D component_histogram(name, name, 
            x_variable_->getBins(), min_x_bin_edge, max_x_bin_edge,
            y_variable_->getBins(), min_y_bin_edge, max_y_bin_edge, 
            z_variable_->getBins(), min_z_bin_edge, max_z_bin_edge);
        component_histograms.push_back(component_histogram);
        name.Clear();
      }
      event_sign_histograms.push_back(component_histograms);
    }
    histograms_.push_back(event_sign_histograms);
  }
}


Fit3D::~Fit3D()
{
  // std::cout << "Starting destructor." << std::endl;
  // TODO member pointers need to be deleted, but these lines give segfault.
  // delete x_variable_;
  // delete y_variable_;
  // delete y_variable_;
}

void Fit3D::ntupleLoopCore(const int& entry_id)
{
  x_value_ = l0_pcm + l1_pcm;
  y_value_ = cos_thta;
  if (evt_sign != 0) {
    z_value_ = abs(l0_ip_dz - l1_ip_dz);
  } else {
    z_value_ = l0_ip_dz - l1_ip_dz;
  }
  
  AnalysisSelectors cuts(*this);
  fillDataSet(cuts.signal_type());
  fillHistograms(cuts.signal_type());  
}

void Fit3D::fillDataSet(const int &component)
{ 
  // WARNING - Missing columns will silently not be added to the data set!
  if (flag_output_a_dataset_ && data_set_) {
    x_variable_->setVal(x_value_);
    y_variable_->setVal(y_value_);
    z_variable_->setVal(z_value_);
    component_->setIndex(component);
    event_sign_->setIndex(evt_sign);
    event_species_->setIndex(typ_tru);
    data_set_->add(
        RooArgSet(
            *x_variable_,
            *y_variable_,
            *z_variable_,
            *component_,
            *event_sign_,
            *event_species_));
  }
}

void Fit3D::fillHistograms(const int& component)
{
  // Catch bad entries. Real data will have typ_true = 0.
  histograms_[typ_tru][evt_sign + 1][component].Fill(
      x_value_,
      y_value_,
      z_value_);
}

void Fit3D::drawHistograms()
{
  for (int i = 0; i < tags_.species_.size(); ++i) {
    TCanvas canvas("canvas", tags_.species_[i], 1200, 1200);
    canvas.Divide(3, 3);
    vector<TH1*> x_histograms;
    vector<TH1*> y_histograms;
    vector<TH1*> z_histograms;
    for (int j = 0; j < tags_.signs_.size(); ++j) {
      x_histograms.clear();
      y_histograms.clear();
      z_histograms.clear();
      int x_max(0);
      int y_max(0);
      int z_max(0);
      for (int k = 0; k < tags_.components_.size(); ++k) {
        TString base_name = histograms_[i][j][k].GetName();
        TH1* x_histogram = histograms_[i][j][k].Project3D("x");
        x_histograms.push_back(x_histogram);
        if (x_histogram->GetMaximum() > x_max) {
          x_max = x_histogram->GetMaximum();
        }
        TH1* y_histogram = histograms_[i][j][k].Project3D("y");
        y_histograms.push_back(y_histogram);
        if (y_histogram->GetMaximum() > y_max) {
          y_max = y_histogram->GetMaximum();
        }
        TH1* z_histogram = histograms_[i][j][k].Project3D("z");
        z_histograms.push_back(z_histogram);
        if (z_histogram->GetMaximum() > z_max) {
          z_max = z_histogram->GetMaximum();
        }
      }
      for (int k = 0; k < tags_.components_.size(); ++k) {
        TString options = ((k == 0) ? "e" : "e same");
        
        canvas.cd(j + 1 + 0);
        x_histograms[k]->SetMaximum(1.1 * x_max);
        x_histograms[k]->SetLineColor(tags_.colors_[k]);
        x_histograms[k]->Draw(options);
        canvas.cd(j + 1 + 3);
        y_histograms[k]->SetMaximum(1.1 * y_max);
        y_histograms[k]->SetLineColor(tags_.colors_[k]);
        y_histograms[k]->Draw(options);
        canvas.cd(j + 1 + 6);
        z_histograms[k]->SetMaximum(1.1 * z_max);
        z_histograms[k]->SetLineColor(tags_.colors_[k]);
        z_histograms[k]->Draw(options);
      }
    }
    canvas.Print(tags_.species_[i] + ".eps");
    while (!x_histograms.empty()) {
      delete x_histograms.back();
      x_histograms.pop_back();
    }
    while (!y_histograms.empty()) {
      delete y_histograms.back();
      y_histograms.pop_back();
    }
    while (!z_histograms.empty()) {
      delete z_histograms.back();
      z_histograms.pop_back();
    }
  }
}

void Fit3D::saveHistograms(const TString& filename)
{
  TFile histogram_file(filename, "RECREATE");
  for (int i = 0; i < tags_.species_.size(); ++i) {
    for (int j = 0; j < tags_.signs_.size(); ++j) {
      for (int k = 0; k < tags_.components_.size(); ++k) {
        histograms_[i][j][k].Write();
      }
    }
  }
  histogram_file.Close();
}

void Fit3D::generateModels(const int& interpolation_order)
{
  std::cout << "Generating models..." << endl;
  vector<TCut> component_cuts;
  component_cuts.push_back(bs_events_cut_);
  component_cuts.push_back(bd_events_cut_);
  component_cuts.push_back(cw_events_cut_);
  component_cuts.push_back(ww_events_cut_);
  component_cuts.push_back(cn_events_cut_);
  vector<TCut> sign_cuts;
  sign_cuts.push_back(nn_events_cut_);
  sign_cuts.push_back(pn_events_cut_);
  sign_cuts.push_back(pp_events_cut_);
  RooWorkspace model_space("model_space", "Fit models");
  RooArgList xyz_variables(*x_variable_, *y_variable_, *z_variable_);
  data_set_->Print();
  for (int j = 0; j < tags_.signs_.size(); j++) {
    for (int k = 1; k < tags_.components_.size(); k++) {
      TString name = "";
      name.Append(tags_.signs_[j]);
      name.Append("_");
      name.Append(tags_.components_[k]);
      std::cout << "Generating " << name << " model." << std::endl;
      TCut cut = sign_cuts[j] + component_cuts[k-1];
      std::cout << cut << std::endl;
      RooDataSet* unbinned_data = (RooDataSet*) data_set_->reduce(
          RooFit::Cut(cut));
      RooDataHist binned_data(name + "_binned_data", name + "_binned_data",
          xyz_variables, *unbinned_data);
      RooHistPdf pdf(name + "_pdf", name + "_pdf",
          xyz_variables, binned_data, interpolation_order);
      model_space.import(pdf);
      if (false) {
        RooPlot* frame = x_variable_->frame();
        unbinned_data->plotOn(frame);
        binned_data.plotOn(frame);
        pdf.plotOn(frame);
        TCanvas c1("c1", "Model", 200, 10, 700, 500);
        frame->Draw();
        c1.Print(name + "_pdf.eps");
        delete frame;
      }
    }
  }
  
  model_space.import(*data_set_);
  TFile models_file("models.root", "RECREATE");
  model_space.Write();
  models_file.Close();
  
  std::cout << "Finished!" << endl;
}

void Fit3D::fitData(const TString& filename, const TString& data_set)
{  
  TFile models_file("models.root", "READ");
  RooWorkspace* model_space = (RooWorkspace*) models_file.Get("model_space");
  
  if (!model_space) {
    std::cout << "Cannot load models. Exiting." << std::endl;
    return;
  }
  
  RooAbsData &fit_data = *data_set_;
  RooAbsPdf* pp_bs_pdf = model_space->pdf("pp_bs_pdf");
  RooAbsPdf* pp_bd_pdf = model_space->pdf("pp_bd_pdf");
  RooAbsPdf* pp_cw_pdf = model_space->pdf("pp_cw_pdf");
  RooAbsPdf* pp_ww_pdf = model_space->pdf("pp_ww_pdf");
  RooAbsPdf* pp_cn_pdf = model_space->pdf("pp_cn_pdf");
  
  RooAbsPdf* nn_bs_pdf = model_space->pdf("nn_bs_pdf");
  RooAbsPdf* nn_bd_pdf = model_space->pdf("nn_bd_pdf");
  RooAbsPdf* nn_cw_pdf = model_space->pdf("nn_cw_pdf");
  RooAbsPdf* nn_ww_pdf = model_space->pdf("nn_ww_pdf");
  RooAbsPdf* nn_cn_pdf = model_space->pdf("nn_cn_pdf");
  
  std::cout << "Generating data sets for fit." << std::endl;
  RooDataSet& pp_data = *((RooDataSet*) fit_data.reduce(pp_events_cut_));
  RooDataSet& nn_data = *((RooDataSet*) fit_data.reduce(nn_events_cut_));
  
  RooRealVar n_bs_pp("n_bs_pp", "n_bs_pp", 0.0000e+00, 1.0000e+6);
  RooRealVar n_bd_pp("n_bd_pp", "n_bd_pp", 0.0000e+00, 1.0000e+6);
  RooRealVar n_cw_pp("n_cw_pp", "n_cw_pp", 0.0000e+00, 1.0000e+6);
  RooRealVar n_ww_pp("n_ww_pp", "n_ww_pp", 0.0000e+00, 1.0000e+6);
  RooRealVar n_cn_pp("n_cn_pp", "n_cn_pp", 0.0000e+00, 1.0000e+6);
  
  RooRealVar n_bs_nn("n_bs_nn", "n_bs_nn", 0.0000e+00, 1.0000e+6);
  RooRealVar n_bd_nn("n_bd_nn", "n_bd_nn", 0.0000e+00, 1.0000e+6);
  RooRealVar n_cw_nn("n_cw_nn", "n_cw_nn", 0.0000e+00, 1.0000e+6);
  RooRealVar n_ww_nn("n_ww_nn", "n_ww_nn", 0.0000e+00, 1.0000e+6);
  RooRealVar n_cn_nn("n_cn_nn", "n_cn_nn", 0.0000e+00, 1.0000e+6);
  
  RooArgList pp_yields(n_bs_pp, n_bd_pp, n_cw_pp, n_ww_pp, n_cn_pp);
  RooArgList nn_yields(n_bs_nn, n_bd_nn, n_cw_nn, n_ww_nn, n_cn_nn);
  
  RooArgList pp_model_components(*pp_bs_pdf, *pp_bd_pdf, *pp_cw_pdf,
      *pp_ww_pdf, *pp_cn_pdf);
  RooArgList nn_model_components(*nn_bs_pdf, *nn_bd_pdf, *nn_cw_pdf,
      *nn_ww_pdf, *nn_cn_pdf);
  
  RooAddPdf pp_model(
      "pp_model", "sig+bak", pp_model_components, pp_yields);
  RooAddPdf nn_model(
      "nn_model", "sig+bak", nn_model_components, nn_yields);
  
  std::cout << "Starting the fit." << std::endl;
  RooFitResult* pp_fit_results = pp_model.fitTo(
      pp_data,
      RooFit::Minos(true),
      RooFit::NumCPU(3),
      RooFit::Timer(true),
      RooFit::Save(true));
  
  RooFitResult* nn_fit_results = nn_model.fitTo(
      nn_data,
      RooFit::Minos(true),
      RooFit::NumCPU(3),
      RooFit::Timer(true),
      RooFit::Save(true));
  
  plotFitAccuracy(pp_data, *pp_fit_results);
  plotFitAccuracy(nn_data, *nn_fit_results);
  plotFitProjection(*x_variable_, pp_data, *pp_fit_results, pp_model,
      *pp_bs_pdf, *pp_bd_pdf, *pp_cw_pdf, *pp_ww_pdf, *pp_cn_pdf,
      "pp_x_fit.eps");
  plotFitProjection(*y_variable_, pp_data, *pp_fit_results, pp_model,
      *pp_bs_pdf, *pp_bd_pdf, *pp_cw_pdf, *pp_ww_pdf, *pp_cn_pdf,
      "pp_y_fit.eps");
  plotFitProjection(*z_variable_, pp_data, *pp_fit_results, pp_model,
      *pp_bs_pdf, *pp_bd_pdf, *pp_cw_pdf, *pp_ww_pdf, *pp_cn_pdf,
      "pp_z_fit.eps");
  plotFitProjection(*x_variable_, nn_data, *nn_fit_results, nn_model,
      *nn_bs_pdf, *nn_bd_pdf, *nn_cw_pdf, *nn_ww_pdf, *nn_cn_pdf,
      "nn_x_fit.eps");
  plotFitProjection(*y_variable_, nn_data, *nn_fit_results, nn_model,
      *nn_bs_pdf, *nn_bd_pdf, *nn_cw_pdf, *nn_ww_pdf, *nn_cn_pdf,
      "nn_y_fit.eps");
  plotFitProjection(*z_variable_, nn_data, *nn_fit_results, nn_model,
      *nn_bs_pdf, *nn_bd_pdf, *nn_cw_pdf, *nn_ww_pdf, *nn_cn_pdf,
      "nn_z_fit.eps");
}

void Fit3D::plotFitAccuracy(
    const RooDataSet& mc_data,
    const RooFitResult& fit)
{
  fit.Print("v");

  double n_true_bs = mc_data.sumEntries("component == component::bs");
  double n_true_bd = mc_data.sumEntries("component == component::bd");
  double n_true_cw = mc_data.sumEntries("component == component::cw");
  double n_true_ww = mc_data.sumEntries("component == component::ww");
  double n_true_cn = mc_data.sumEntries("component == component::cn");

  RooRealVar* bs_fit = (RooRealVar*) fit.floatParsFinal().find("n_bs_pp");
  RooRealVar* bd_fit = (RooRealVar*) fit.floatParsFinal().find("n_bd_pp");
  RooRealVar* cw_fit = (RooRealVar*) fit.floatParsFinal().find("n_cw_pp");
  RooRealVar* ww_fit = (RooRealVar*) fit.floatParsFinal().find("n_ww_pp");
  RooRealVar* cn_fit = (RooRealVar*) fit.floatParsFinal().find("n_cn_pp");
  
  TString title("Fit Accuracy (N^{++}_{fit}-N^{++}_{true})");
  TString filename("fit_accuracy_pp.eps");
  if (!bs_fit) {
    bs_fit = (RooRealVar*) fit.floatParsFinal().find("n_bs_nn");
    bd_fit = (RooRealVar*) fit.floatParsFinal().find("n_bd_nn");
    cw_fit = (RooRealVar*) fit.floatParsFinal().find("n_cw_nn");
    ww_fit = (RooRealVar*) fit.floatParsFinal().find("n_ww_nn");
    cn_fit = (RooRealVar*) fit.floatParsFinal().find("n_cn_nn");
    title = TString("Fit Accuracy (N^{--}_{fit}-N^{--}_{true})");
    filename = TString("fit_accuracy_nn.eps");
  }
  if (!bs_fit) {
    // Error. Quit while ahead.
    cout << "Error in plotFitAccuracy(): "
         << "Cannot find fit variables. Check names are valid."
         << endl;
    return;
  }

  std::cout << n_true_bd << std::endl;
  std::cout << n_true_bs << std::endl;
  std::cout << n_true_cn << std::endl;
  std::cout << n_true_cw << std::endl;
  std::cout << n_true_ww << std::endl;
  
  TCanvas c1("c1", title, 200, 10, 700, 500);
  c1.SetGrid();
  double x[5] = {1, 2, 3, 4, 5};
  double y[5] = {
      bs_fit->getVal() - n_true_bs,
      bd_fit->getVal() - n_true_bd,
      cw_fit->getVal() - n_true_cw,
      ww_fit->getVal() - n_true_ww,
      cn_fit->getVal() - n_true_cn};
  double exl[5] = {0, 0, 0, 0, 0};
  double exh[5] = {0, 0, 0, 0, 0};
  double eyl[5] = {
      -bs_fit->getErrorLo(),
      -bd_fit->getErrorLo(),
      -cw_fit->getErrorLo(),
      -ww_fit->getErrorLo(),
      -cn_fit->getErrorLo()};
  double eyh[5] = {
      bs_fit->getErrorHi(),
      bd_fit->getErrorHi(),
      cw_fit->getErrorHi(),
      ww_fit->getErrorHi(),
      cn_fit->getErrorHi()};
  TGraphAsymmErrors* gr = new TGraphAsymmErrors(5, x, y, exl, exh, eyl, eyh);
  
  TLatex* cc_bs_label = new TLatex(gr->GetX()[0], gr->GetY()[0], " CC (B_{s})");
  TLatex* cc_bd_label = new TLatex(gr->GetX()[1], gr->GetY()[1], " CC (B_{d})");
  TLatex* cw_label = new TLatex(gr->GetX()[2], gr->GetY()[2], " CW");
  TLatex* ww_label = new TLatex(gr->GetX()[3], gr->GetY()[3], " WW");
  TLatex* cn_label = new TLatex(gr->GetX()[4], gr->GetY()[4], " CN");
  
  gr->GetListOfFunctions()->Add(cc_bs_label);
  gr->GetListOfFunctions()->Add(cc_bd_label);
  gr->GetListOfFunctions()->Add(cw_label);
  gr->GetListOfFunctions()->Add(ww_label);
  gr->GetListOfFunctions()->Add(cn_label);
  
  gr->SetTitle(title);
  gr->SetMarkerStyle(kOpenCircle);
  gr->SetMarkerColor(4);
  gr->Draw("AP");

  c1.Print(filename);
  return;
}

void Fit3D::plotFitProjection(
    const RooRealVar &independant_variable,
    const RooDataSet &data,
    const RooFitResult& fit,
    const RooAbsPdf &model,
    const RooAbsPdf &bs_pdf,
    const RooAbsPdf &bd_pdf,
    const RooAbsPdf &cw_pdf,
    const RooAbsPdf &ww_pdf,
    const RooAbsPdf &cn_pdf,
    const TString &filename)
{
  RooPlot* frame = independant_variable.frame();
  TString frame_title = "Fit Projection on ";
  frame_title.Append(independant_variable.GetTitle());
  frame->SetTitle(frame_title);
  
  RooRealVar* bs_fit = (RooRealVar*) fit.floatParsFinal().find("n_bs_pp");
  RooRealVar* bd_fit = (RooRealVar*) fit.floatParsFinal().find("n_bd_pp");
  RooRealVar* cw_fit = (RooRealVar*) fit.floatParsFinal().find("n_cw_pp");
  RooRealVar* ww_fit = (RooRealVar*) fit.floatParsFinal().find("n_ww_pp");
  RooRealVar* cn_fit = (RooRealVar*) fit.floatParsFinal().find("n_cn_pp");
  
  if (!bs_fit) {
    bs_fit = (RooRealVar*) fit.floatParsFinal().find("n_bs_nn");
    bd_fit = (RooRealVar*) fit.floatParsFinal().find("n_bd_nn");
    cw_fit = (RooRealVar*) fit.floatParsFinal().find("n_cw_nn");
    ww_fit = (RooRealVar*) fit.floatParsFinal().find("n_ww_nn");
    cn_fit = (RooRealVar*) fit.floatParsFinal().find("n_cn_nn");
  }
  if (!bs_fit) {
    // Error. Quit while ahead.
    cout << "Error in plotFitAccuracy(): "
         << "Cannot find fit variables. Check names are valid."
         << endl;
    return;
  }
  
  data.plotOn(frame, RooFit::Name("data"));
  // model.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kBlue));
  bs_pdf.plotOn(frame,
      RooFit::Normalization(bs_fit->getVal(), RooAbsReal::NumEvent),
      RooFit::LineStyle(kDashed),
      RooFit::LineWidth(1),
      RooFit::LineColor(kYellow + 2));
  data.plotOn(frame, RooFit::Cut(bs_events_cut_),
      RooFit::LineColor(kYellow),
      RooFit::MarkerStyle(kFullDotMedium));
  
  bd_pdf.plotOn(frame,
      RooFit::Normalization(bd_fit->getVal(), RooAbsReal::NumEvent),
      RooFit::LineStyle(kDashed),
      RooFit::LineWidth(1),
      RooFit::LineColor(kRed + 2));
  data.plotOn(frame, RooFit::Cut(bd_events_cut_),
      RooFit::LineColor(kRed),
      RooFit::MarkerStyle(kFullDotMedium));
  cw_pdf.plotOn(frame,
      RooFit::Normalization(cw_fit->getVal(), RooAbsReal::NumEvent),
      RooFit::LineStyle(kDashed),
      RooFit::LineWidth(1),
      RooFit::LineColor(kGreen + 2));
  data.plotOn(frame, RooFit::Cut(cw_events_cut_),
      RooFit::LineColor(kGreen),
      RooFit::MarkerStyle(kFullDotMedium));
  ww_pdf.plotOn(frame,
      RooFit::Normalization(ww_fit->getVal(), RooAbsReal::NumEvent),
      RooFit::LineStyle(kDashed),
      RooFit::LineWidth(1),
      RooFit::LineColor(kBlue + 2));
  data.plotOn(frame, RooFit::Cut(ww_events_cut_),
      RooFit::LineColor(kBlue),
      RooFit::MarkerStyle(kFullDotMedium));
  cn_pdf.plotOn(frame,
      RooFit::Normalization(cn_fit->getVal(), RooAbsReal::NumEvent),
      RooFit::LineStyle(kDashed),
      RooFit::LineWidth(1),
      RooFit::LineColor(kCyan + 2));
  data.plotOn(frame, RooFit::Cut(cn_events_cut_),
      RooFit::LineColor(kCyan),
      RooFit::MarkerStyle(kFullDotMedium));
      
  TCanvas* c1 = new TCanvas("c1", "Projection", 200, 10, 700, 500);
  frame->Draw();
  c1->Print(filename);
}