#include "Fit3D.h"


void plotFitProjection(const RooRealVar &independant_variable, const RooDataSet &data,
        const RooAbsPdf &model, const TString &filename) {
  RooPlot* frame = independant_variable.frame();
  data.plotOn(frame, RooFit::Name("data"));
  model.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kBlue));
  model.plotOn(frame, RooFit::Components("*bs*"),
      RooFit::LineStyle(kDashed), RooFit::LineWidth(1), RooFit::LineColor(kYellow));
  model.plotOn(frame, RooFit::Components("*bd*"),
      RooFit::LineStyle(kDashed), RooFit::LineWidth(1), RooFit::LineColor(kRed));
  model.plotOn(frame, RooFit::Components("*cw*"),
      RooFit::LineStyle(kDashed), RooFit::LineWidth(1), RooFit::LineColor(kGreen));
  model.plotOn(frame, RooFit::Components("*ww*"),
      RooFit::LineStyle(kDashed), RooFit::LineWidth(1), RooFit::LineColor(kBlue));
  model.plotOn(frame, RooFit::Components("*cn*"),
      RooFit::LineStyle(kDashed), RooFit::LineWidth(1), RooFit::LineColor(kCyan));

  TCanvas* c1 = new TCanvas("c1", "Projection", 200, 10, 700, 500);
  frame->Draw();
  c1->Print(filename);
}

Fit3D::Fit3D(
    TString input_ntuple_file,
    TString analysis_name,
    TString x_axis_label,
    double min_x_bin_edge,
    double max_x_bin_edge,
    TString y_axis_label,
    double min_y_bin_edge,
    double max_y_bin_edge,
    TString z_axis_label,
    double min_z_bin_edge,
    double max_z_bin_edge)
  : DileptonEvents(input_ntuple_file, analysis_name)
{
  std::cout << "Initiliazing Fit3D class." << std::endl;
  
  // Define the data set elements and add them to the dataset.
  x_variable_ = new RooRealVar(
      "x_variable", x_axis_label, min_x_bin_edge, max_x_bin_edge);
  y_variable_ = new RooRealVar(
      "y_variable", y_axis_label, min_y_bin_edge, max_y_bin_edge);
  z_variable_ = new RooRealVar(
      "z_variable", z_axis_label, min_z_bin_edge, max_z_bin_edge);
  x_variable_->setBins(100);
  y_variable_->setBins(50);
  z_variable_->setBins(50);
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
            100, min_x_bin_edge, max_x_bin_edge,
            100, min_y_bin_edge, max_y_bin_edge, 
            100, min_z_bin_edge, max_z_bin_edge);
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

void Fit3D::ntupleLoopCore()
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
    vector<TH1D*> x_histograms;
    vector<TH1D*> y_histograms;
    vector<TH1D*> z_histograms;
    for (int j = 0; j < tags_.signs_.size(); ++j) {
      x_histograms.clear();
      y_histograms.clear();
      z_histograms.clear();
      int x_max(0);
      int y_max(0);
      int z_max(0);
      for (int k = 0; k < tags_.components_.size(); ++k) {
        TString base_name = histograms_[i][j][k].GetName();
        TH1D* x_histogram = histograms_[i][j][k].ProjectionX(base_name + "_x");
        x_histograms.push_back(x_histogram);
        if (x_histogram->GetMaximum() > x_max) {
          x_max = x_histogram->GetMaximum();
        }
        TH1D* y_histogram = histograms_[i][j][k].ProjectionY(base_name + "_y");
        y_histograms.push_back(y_histogram);
        if (y_histogram->GetMaximum() > y_max) {
          y_max = y_histogram->GetMaximum();
        }
        TH1D* z_histogram = histograms_[i][j][k].ProjectionZ(base_name + "_z");
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

void Fit3D::generateModels()
{
  std::cout << "Generating models..." << endl;

  // Reduce data to components.
  RooDataSet* bs_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(bs_events_cut_ && pp_events_cut_));
  RooDataSet* bd_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(bd_events_cut_ && pp_events_cut_));
  RooDataSet* cw_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cw_events_cut_ && pp_events_cut_));
  RooDataSet* ww_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(ww_events_cut_ && pp_events_cut_));
  RooDataSet* cn_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cn_events_cut_ && pp_events_cut_));
  
  /*
  RooDataSet* bs_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(bs_events_cut_ && pn_events_cut_));
  RooDataSet* bd_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(bd_events_cut_ && pn_events_cut_));
  RooDataSet* cw_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cw_events_cut_ && pn_events_cut_));
  RooDataSet* ww_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(ww_events_cut_ && pn_events_cut_));
  RooDataSet* cn_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cn_events_cut_ && pn_events_cut_));
  */
  
  RooDataSet* bs_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(bs_events_cut_ && nn_events_cut_));
  RooDataSet* bd_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(bd_events_cut_ && nn_events_cut_));
  RooDataSet* cw_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cw_events_cut_ && nn_events_cut_));
  RooDataSet* ww_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(ww_events_cut_ && nn_events_cut_));
  RooDataSet* cn_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cn_events_cut_ && nn_events_cut_));

  RooArgList observables(*x_variable_, *y_variable_, *z_variable_);
  RooDataHist bs_pp_binned_data("bs_pp_binned_data", "bs_pp_binned_data",
      observables, *bs_pp_data);
  RooDataHist bd_pp_binned_data("bd_pp_binned_data", "bd_pp_binned_data",
      observables, *bd_pp_data);
  RooDataHist cw_pp_binned_data("cw_pp_binned_data", "cw_pp_binned_data",
      observables, *cw_pp_data);
  RooDataHist ww_pp_binned_data("ww_pp_binned_data", "ww_pp_binned_data",
      observables, *ww_pp_data);
  RooDataHist cn_pp_binned_data("cn_pp_binned_data", "cn_pp_binned_data",
      observables, *cn_pp_data);
  RooDataHist bs_nn_binned_data("bs_nn_binned_data", "bs_nn_binned_data",
      observables, *bs_nn_data);
  RooDataHist bd_nn_binned_data("bd_nn_binned_data", "bd_nn_binned_data",
      observables, *bd_nn_data);
  RooDataHist cw_nn_binned_data("cw_nn_binned_data", "cw_nn_binned_data",
      observables, *cw_nn_data);
  RooDataHist ww_nn_binned_data("ww_nn_binned_data", "ww_nn_binned_data",
      observables, *ww_nn_data);
  RooDataHist cn_nn_binned_data("cn_nn_binned_data", "cn_nn_binned_data",
      observables, *cn_nn_data);

  int interpolation_order = 0;
  RooHistPdf bs_pp_pdf("bs_pp_pdf", "bs_pp_pdf",
      observables, bs_pp_binned_data, interpolation_order);
  RooHistPdf bd_pp_pdf("bd_pp_pdf", "bd_pp_pdf",
      observables, bd_pp_binned_data, interpolation_order);
  RooHistPdf cw_pp_pdf("cw_pp_pdf", "cw_pp_pdf",
      observables, cw_pp_binned_data, interpolation_order);
  RooHistPdf ww_pp_pdf("ww_pp_pdf", "ww_pp_pdf",
      observables, ww_pp_binned_data, interpolation_order);
  RooHistPdf cn_pp_pdf("cn_pp_pdf", "cn_pp_pdf",
      observables, cn_pp_binned_data, interpolation_order);
  RooHistPdf bs_nn_pdf("bs_nn_pdf", "bs_nn_pdf",
      observables, bs_nn_binned_data, interpolation_order);
  RooHistPdf bd_nn_pdf("bd_nn_pdf", "bd_nn_pdf",
      observables, bd_nn_binned_data, interpolation_order);
  RooHistPdf cw_nn_pdf("cw_nn_pdf", "cw_nn_pdf",
      observables, cw_nn_binned_data, interpolation_order);
  RooHistPdf ww_nn_pdf("ww_nn_pdf", "ww_nn_pdf",
      observables, ww_nn_binned_data, interpolation_order);
  RooHistPdf cn_nn_pdf("cn_nn_pdf", "cn_nn_pdf",
      observables, cn_nn_binned_data, interpolation_order);
  
  RooWorkspace model_space("model_space", "Fit models");
  model_space.import(*data_set_);
  model_space.import(bs_pp_pdf);
  model_space.import(bd_pp_pdf);
  model_space.import(cw_pp_pdf);
  model_space.import(ww_pp_pdf);
  model_space.import(cn_pp_pdf);
  model_space.import(bs_nn_pdf);
  model_space.import(bd_nn_pdf);
  model_space.import(cw_nn_pdf);
  model_space.import(ww_nn_pdf);
  model_space.import(cn_nn_pdf);
  TFile models_file("models.root", "RECREATE");
  model_space.Write();
  models_file.Close();
  
  std::cout << "Finished!" << endl;
}

void Fit3D::fitData(const TString& filename, const TString& data_set)
{

  RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  
  TFile models_file("models.root", "READ");
  RooWorkspace* model_space = (RooWorkspace*) models_file.Get("model_space");
  
  if (!model_space) {
    std::cout << "Cannot load models. Exiting." << std::endl;
    return;
  }
  
  RooAbsData &all_data = *(model_space->data("data_set"));
  RooAbsPdf* bs_pp_pdf = model_space->pdf("bs_pp_pdf");
  RooAbsPdf* bd_pp_pdf = model_space->pdf("bd_pp_pdf");
  RooAbsPdf* cw_pp_pdf = model_space->pdf("cw_pp_pdf");
  RooAbsPdf* ww_pp_pdf = model_space->pdf("ww_pp_pdf");
  RooAbsPdf* cn_pp_pdf = model_space->pdf("cn_pp_pdf");
  RooAbsPdf* bs_nn_pdf = model_space->pdf("bs_nn_pdf");
  RooAbsPdf* bd_nn_pdf = model_space->pdf("bd_nn_pdf");
  RooAbsPdf* cw_nn_pdf = model_space->pdf("cw_nn_pdf");
  RooAbsPdf* ww_nn_pdf = model_space->pdf("ww_nn_pdf");
  RooAbsPdf* cn_nn_pdf = model_space->pdf("cn_nn_pdf");
  
  std::cout << "Generating data sets for fit." << std::endl;
  RooDataSet& pp_data = *((RooDataSet*) all_data.reduce(pp_events_cut_));
  RooDataSet& nn_data = *((RooDataSet*) all_data.reduce(nn_events_cut_));
  
  RooRealVar n_bs_pp("n_bs_pp", "n_bs_pp", 0.0000e+00, 1.0000e+10);
  RooRealVar n_bd_pp("n_bd_pp", "n_bd_pp", 0.0000e+00, 1.0000e+10);
  RooRealVar n_cw_pp("n_cw_pp", "n_cw_pp", 0.0000e+00, 1.0000e+10);
  RooRealVar n_ww_pp("n_ww_pp", "n_ww_pp", 0.0000e+00, 1.0000e+10);
  RooRealVar n_cn_pp("n_cn_pp", "n_cn_pp", 0.0000e+00, 1.0000e+10);
  
  /*
  RooRealVar n_bs_pn("n_bs_pn", "n_bs_pn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_bd_pn("n_bd_pn", "n_bd_pn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_cw_pn("n_cw_pn", "n_cw_pn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_ww_pn("n_ww_pn", "n_ww_pn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_cn_pn("n_cn_pn", "n_cn_pn", 0.0000e+00, 1.0000e+10);
  */
  
  RooRealVar n_bs_nn("n_bs_nn", "n_bs_nn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_bd_nn("n_bd_nn", "n_bd_nn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_cw_nn("n_cw_nn", "n_cw_nn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_ww_nn("n_ww_nn", "n_ww_nn", 0.0000e+00, 1.0000e+10);
  RooRealVar n_cn_nn("n_cn_nn", "n_cn_nn", 0.0000e+00, 1.0000e+10);
  
  RooArgList pp_yields(n_bs_pp, n_bd_pp, n_cw_pp, n_ww_pp, n_cn_pp);
  /*
  RooArgList pn_yields(n_bs_pn, n_bd_pn, n_cw_pn, n_ww_pn, n_cn_pn);
  */
  RooArgList nn_yields(n_bs_nn, n_bd_nn, n_cw_nn, n_ww_nn, n_cn_nn);
  
  RooArgList pp_model_components(
      *bs_pp_pdf, *bd_pp_pdf, *cw_pp_pdf, *ww_pp_pdf, *cn_pp_pdf);
  RooArgList nn_model_components(
      *bs_nn_pdf, *bd_nn_pdf, *cw_nn_pdf, *ww_nn_pdf, *cn_nn_pdf);
  
  RooAddPdf pp_model("pp_model", "sig+bak", pp_model_components, pp_yields);
  RooAddPdf nn_model("nn_model", "sig+bak", nn_model_components, nn_yields);
  
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
  plotFitProjection(*x_variable_, pp_data, pp_model, "pp_x_fit.eps");
  plotFitProjection(*y_variable_, pp_data, pp_model, "pp_y_fit.eps");
  plotFitProjection(*z_variable_, pp_data, pp_model, "pp_z_fit.eps");
  plotFitProjection(*x_variable_, nn_data, nn_model, "nn_x_fit.eps");
  plotFitProjection(*y_variable_, nn_data, nn_model, "nn_y_fit.eps");
  plotFitProjection(*z_variable_, nn_data, nn_model, "nn_z_fit.eps");
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
    cout << "Error in plotFitAccuracy(): Cannot find fit variables. Check names are valid."
         << endl;
    return;
  }

  TCanvas* c1 = new TCanvas("c1", title, 200, 10, 700, 500);
  c1->SetGrid();
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
  gr->GetXaxis()->SetBinLabel(1, "CC B_{s}");
  gr->GetXaxis()->SetBinLabel(1, "CC B_{d}");
  gr->GetXaxis()->SetBinLabel(1, "CW");
  gr->GetXaxis()->SetBinLabel(1, "WW");
  gr->GetXaxis()->SetBinLabel(1, "CN");
  gr->SetTitle(title);
  gr->SetMarkerStyle(kOpenCircle);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");

  c1->Print(filename);
  return;
}