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
  
  output_path_
      = "/home/matt/research/belle/adcab/analysis/DileptonEvents/Fit3D/output/"
      + analysis_name + "/";
  
  boost::filesystem::path output_dir(output_path_);
  
  if (boost::filesystem::exists(output_dir)) {
    cout << "Output path exists. Overwriting directory contents!\n";
  } else {
    cout << "Output path does not exist. Creating new directory.\n";
    boost::filesystem::create_directories(output_dir);
  }
  
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
  x_variable_->setBins(40);
  y_variable_->setBins(40);
  z_variable_->setBins(40);
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
    event_species_->setIndex(typ_asn);
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
  // Use typ_asn for real data.
  histograms_[typ_asn][evt_sign + 1][component].Fill(
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
      int x_min(0);
      int y_min(0);
      int z_min(0);
      for (int k = 0; k < tags_.components_.size(); ++k) {
        TString base_name = histograms_[i][j][k].GetName();
        TH1* x_histogram = histograms_[i][j][k].Project3D("x");
        x_histograms.push_back(x_histogram);
        if (x_histogram->GetMinimum() > x_min) {
          x_min = x_histogram->GetMinimum();
        }
        TH1* y_histogram = histograms_[i][j][k].Project3D("y");
        y_histograms.push_back(y_histogram);
        if (y_histogram->GetMinimum() > y_min) {
          y_min = y_histogram->GetMinimum();
        }
        TH1* z_histogram = histograms_[i][j][k].Project3D("z");
        z_histograms.push_back(z_histogram);
        if (z_histogram->GetMinimum() > z_min) {
          z_min = z_histogram->GetMinimum();
        }
      }
      
      TH1F* x_total = (TH1F*) x_histograms[0]->Clone("X Total");
      TH1F* y_total = (TH1F*) y_histograms[0]->Clone("Y Total");
      TH1F* z_total = (TH1F*) z_histograms[0]->Clone("Z Total");
      
      for (int k = 1; k < tags_.components_.size(); ++k) {
        x_total->Add(x_histograms[k]);
        y_total->Add(y_histograms[k]);
        z_total->Add(z_histograms[k]);
      }
    
      canvas.cd(j + 1 + 0);
      x_total->SetLineColor(1);
      x_total->SetMinimum(0);
      x_total->Draw("e");
      canvas.cd(j + 1 + 3);
      y_total->SetLineColor(1);
      y_total->SetMinimum(0);
      y_total->Draw("e");
      canvas.cd(j + 1 + 6);
      z_total->SetLineColor(1);
      z_total->SetMinimum(0);
      z_total->Draw("e");
      
      for (int k = 0; k < tags_.components_.size(); ++k) {
        TString options = "e same";
        
        canvas.cd(j + 1 + 0);
        x_histograms[k]->SetLineColor(tags_.colors_[k]);
        x_histograms[k]->Draw(options);
        canvas.cd(j + 1 + 3);
        y_histograms[k]->SetLineColor(tags_.colors_[k]);
        y_histograms[k]->Draw(options);
        canvas.cd(j + 1 + 6);
        z_histograms[k]->SetLineColor(tags_.colors_[k]);
        z_histograms[k]->Draw(options);
      }
    }
    canvas.Print(output_path_ + tags_.species_[i] + ".eps");
    
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
  RooArgList xy_variables(*x_variable_, *y_variable_);
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
      RooDataHist binned_xy_data(name + "_binned_xy_data", name + "_binned_xy_data",
          xy_variables, *unbinned_data);
      RooDataHist binned_z_data(name + "_binned_z_data", name + "_binned_z_data",
          *z_variable_, *unbinned_data);
      RooHistPdf xy_pdf(name + "_xy_pdf", name + "_xy_pdf",
          xy_variables, binned_xy_data, interpolation_order);
      RooHistPdf z_pdf(name + "_z_pdf", name + "_z_pdf",
          *z_variable_, binned_z_data, interpolation_order);
      RooProdPdf pdf(name + "_pdf", name + "_pdf", xy_pdf, z_pdf);
      model_space.import(pdf);
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
  
  std::cout << "Organizing data sets for fit." << std::endl;
  RooDataSet& pp_data = *((RooDataSet*) fit_data.reduce(pp_events_cut_));
  RooDataSet& nn_data = *((RooDataSet*) fit_data.reduce(nn_events_cut_));
  
  double pp_data_size = double(pp_data.numEntries());
  double nn_data_size = double(pp_data.numEntries());
  
  RooRealVar n_bs_pp("n_bs_pp", "n_bs_pp", 200, 0, 50 * pp_data_size);
  RooRealVar n_bd_pp("n_bd_pp", "n_bd_pp", 200, 0, 50 * pp_data_size);
  RooRealVar n_cw_pp("n_cw_pp", "n_cw_pp", 200, 0, 50 * pp_data_size);
  RooRealVar n_ww_pp("n_ww_pp", "n_ww_pp", 200, 0, 50 * pp_data_size);
  RooRealVar n_cn_pp("n_cn_pp", "n_cn_pp", 200, 0, 50 * pp_data_size);
  
  RooRealVar n_bs_nn("n_bs_nn", "n_bs_nn", 200, 0, 50 * nn_data_size);
  RooRealVar n_bd_nn("n_bd_nn", "n_bd_nn", 200, 0, 50 * nn_data_size);
  RooRealVar n_cw_nn("n_cw_nn", "n_cw_nn", 200, 0, 50 * nn_data_size);
  RooRealVar n_ww_nn("n_ww_nn", "n_ww_nn", 200, 0, 50 * nn_data_size);
  RooRealVar n_cn_nn("n_cn_nn", "n_cn_nn", 200, 0, 50 * nn_data_size);
  
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
  pp_fit_results_ = pp_model.fitTo(
      pp_data,
      RooFit::Minos(true),
      RooFit::Hesse(false),
      RooFit::NumCPU(3),
      RooFit::Timer(true),
      RooFit::Save(true));
  
  nn_fit_results_ = nn_model.fitTo(
      nn_data,
      RooFit::Minos(true),
      RooFit::Hesse(false),
      RooFit::NumCPU(3),
      RooFit::Timer(true),
      RooFit::Save(true));
  
  plotFitAccuracy(pp_data, *pp_fit_results_);
  plotFitAccuracy(nn_data, *nn_fit_results_);
  plotFitProjection(*x_variable_, pp_data, *pp_fit_results_, pp_model,
      *pp_bs_pdf, *pp_bd_pdf, *pp_cw_pdf, *pp_ww_pdf, *pp_cn_pdf,
      "pp_x_fit.eps");
  plotFitProjection(*y_variable_, pp_data, *pp_fit_results_, pp_model,
      *pp_bs_pdf, *pp_bd_pdf, *pp_cw_pdf, *pp_ww_pdf, *pp_cn_pdf,
      "pp_y_fit.eps");
  plotFitProjection(*z_variable_, pp_data, *pp_fit_results_, pp_model,
      *pp_bs_pdf, *pp_bd_pdf, *pp_cw_pdf, *pp_ww_pdf, *pp_cn_pdf,
      "pp_z_fit.eps");
  plotFitProjection(*x_variable_, nn_data, *nn_fit_results_, nn_model,
      *nn_bs_pdf, *nn_bd_pdf, *nn_cw_pdf, *nn_ww_pdf, *nn_cn_pdf,
      "nn_x_fit.eps");
  plotFitProjection(*y_variable_, nn_data, *nn_fit_results_, nn_model,
      *nn_bs_pdf, *nn_bd_pdf, *nn_cw_pdf, *nn_ww_pdf, *nn_cn_pdf,
      "nn_y_fit.eps");
  plotFitProjection(*z_variable_, nn_data, *nn_fit_results_, nn_model,
      *nn_bs_pdf, *nn_bd_pdf, *nn_cw_pdf, *nn_ww_pdf, *nn_cn_pdf,
      "nn_z_fit.eps");
}

void Fit3D::plotAsymmetry() {
  RooRealVar* n_bs_pp
      = (RooRealVar*) pp_fit_results_->floatParsFinal().find("n_bs_pp");
  RooRealVar* n_bd_pp
      = (RooRealVar*) pp_fit_results_->floatParsFinal().find("n_bd_pp");
  
  RooRealVar* n_bs_nn
      = (RooRealVar*) nn_fit_results_->floatParsFinal().find("n_bs_nn");
  RooRealVar* n_bd_nn
      = (RooRealVar*) nn_fit_results_->floatParsFinal().find("n_bd_nn");
  
  double bs_population_difference(n_bs_pp->getVal() - n_bs_nn->getVal());
  double bd_population_difference(n_bd_pp->getVal() - n_bd_nn->getVal());
  double bs_population_sum(n_bs_pp->getVal() + n_bs_nn->getVal());
  double bd_population_sum(n_bd_pp->getVal() + n_bd_nn->getVal());
  
  double error_n_bs_pp = n_bs_pp->getErrorLo();
  if  (n_bs_pp->getErrorHi() > error_n_bs_pp) {
    error_n_bs_pp = n_bs_pp->getErrorHi();
  }
  double error_n_bs_nn = n_bs_nn->getErrorLo();
  if  (n_bs_nn->getErrorHi() > error_n_bs_nn) {
    error_n_bs_nn = n_bs_nn->getErrorHi();
  }
  
  double error_n_bd_pp = n_bd_pp->getErrorLo();
  if  (n_bd_pp->getErrorHi() > error_n_bd_pp) {
    error_n_bd_pp = n_bd_pp->getErrorHi();
  }
  double error_n_bd_nn = n_bd_nn->getErrorLo();
  if  (n_bd_nn->getErrorHi() > error_n_bd_nn) {
    error_n_bd_nn = n_bd_nn->getErrorHi();
  }
  
  double n_b_pp = n_bs_pp->getVal() + n_bd_pp->getVal();
  double n_b_nn = n_bs_nn->getVal() + n_bd_nn->getVal();
  
  double error_n_b_pp = sqrt(error_n_bs_pp * error_n_bs_pp
      + error_n_bd_pp * error_n_bd_pp);
  double error_n_b_nn = sqrt(error_n_bs_nn * error_n_bs_nn
      + error_n_bd_nn * error_n_bd_nn);
  
  double bs_asymmetry = bs_population_difference / bs_population_sum;
  double bs_asymmetry_error = 2 * sqrt(
      (pow(error_n_bs_pp, 2) * pow(n_bs_nn->getVal(), 2)
      + pow(error_n_bs_nn, 2) * pow(n_bs_pp->getVal(), 2))
      / pow(bs_population_sum, 4));
  
  double bd_asymmetry = bd_population_difference / bd_population_sum;
  double bd_asymmetry_error = 2 * sqrt(
      (pow(error_n_bd_pp, 2) * pow(n_bd_nn->getVal(), 2)
      + pow(error_n_bd_nn, 2) * pow(n_bd_pp->getVal(), 2))
      / pow(bd_population_sum, 4));
  
  double b_asymmetry = (bs_population_difference + bd_population_difference)
      / (bs_population_sum + bd_population_sum);
  double b_asymmetry_error = 2 * sqrt(
      (pow(error_n_b_pp, 2) * pow(n_b_nn, 2)
      + pow(error_n_b_nn, 2) * pow(n_b_pp, 2)) / pow(n_b_pp + n_b_nn, 4));
  
  double asymmetries[3] = {
      bs_asymmetry,
      bd_asymmetry,
      b_asymmetry};
  double asymmetry_errors[3] = {
      bs_asymmetry_error,
      bd_asymmetry_error,
      b_asymmetry_error};
  double x[3] = {1, 2, 3};
  double x_errors[3] = {0, 0, 0};
  
  TGraph line(2);
  line.SetPoint(0,0,0);
  line.SetPoint(1,7,0);
  line.SetLineColor(kRed);
  
  TGraphErrors asymmetries_graph(3, x, asymmetries, x_errors, asymmetry_errors);
  asymmetries_graph.SetMarkerStyle(kOpenCircle);
  asymmetries_graph.SetTitle("Asymmetries");
  asymmetries_graph.GetYaxis()->SetTitle("Asymmetry");
  asymmetries_graph.GetXaxis()->SetTitle("a_{s}, a_{d}, A_{sl}^{b}");
  
  TCanvas c1("Asymmetries", "Asymmetries", 800, 600);
  c1.cd(1);
  asymmetries_graph.Draw("AP");
  line.Draw("l");
  c1.Print(output_path_ + "asymmetries.eps");
  
  TFile plots_file(output_path_ + "asymmetries.root", "RECREATE");
  asymmetries_graph.Write();
  plots_file.Close();
  
  return;
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
  TString filename(output_path_ + "fit_accuracy_pp");
  if (!bs_fit) {
    bs_fit = (RooRealVar*) fit.floatParsFinal().find("n_bs_nn");
    bd_fit = (RooRealVar*) fit.floatParsFinal().find("n_bd_nn");
    cw_fit = (RooRealVar*) fit.floatParsFinal().find("n_cw_nn");
    ww_fit = (RooRealVar*) fit.floatParsFinal().find("n_ww_nn");
    cn_fit = (RooRealVar*) fit.floatParsFinal().find("n_cn_nn");
    title = TString("Fit Accuracy (N^{--}_{fit}-N^{--}_{true})");
    filename = TString(output_path_ + "fit_accuracy_nn");
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

  c1.Print(filename + ".eps");
  
  TFile accuracy_plot_file(filename + ".root", "RECREATE");
  gr->Write();
  accuracy_plot_file.Close();
  
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
  model.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kBlue));
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
  c1->Print(output_path_ + filename);
}