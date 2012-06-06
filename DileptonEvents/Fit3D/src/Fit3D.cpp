#include "Fit3D.h"
#include "EventSelectors.h"
#include <iostream>
#include "TCut.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"


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
        HistogramListOverVariables variable_histograms;
        variable_histograms.clear();
        name.Append(tags_.species_[i]);
        name.Append("_");
        name.Append(tags_.signs_[j]);
        name.Append("_");
        name.Append(tags_.components_[k]);
        TH1F x_histogram(name + "_x", name + "_x",
            100, min_x_bin_edge, max_x_bin_edge);
        TH1F y_histogram(name + "_y", name + "_y",
            100, min_y_bin_edge, max_y_bin_edge);
        TH1F z_histogram(name + "_z", name + "_z",
            100, min_z_bin_edge, max_z_bin_edge);
        variable_histograms.push_back(x_histogram);
        variable_histograms.push_back(y_histogram);
        variable_histograms.push_back(z_histogram);
        component_histograms.push_back(variable_histograms);
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
  histograms_[typ_tru][evt_sign + 1][component][0].Fill(x_value_);
  histograms_[typ_tru][evt_sign + 1][component][1].Fill(y_value_);
  histograms_[typ_tru][evt_sign + 1][component][2].Fill(z_value_);
}

void Fit3D::drawHistograms()
{
  for (int i = 0; i < tags_.species_.size(); ++i) {
    TCanvas canvas("canvas", tags_.species_[i], 1200, 1200);
    canvas.Divide(3, 3);
    int canvas_column(0);
    for (int j = 0; j < tags_.signs_.size(); ++j) {
      double greatest_x_bin_global(0);
      double greatest_y_bin_global(0);
      double greatest_z_bin_global(0);
      for (int k = 0; k < tags_.components_.size(); ++k) {
        double greatest_x_bin_local = histograms_[i][j][k][0].GetMaximum();
        double greatest_y_bin_local = histograms_[i][j][k][1].GetMaximum();
        double greatest_z_bin_local = histograms_[i][j][k][2].GetMaximum();
        if (greatest_x_bin_local > greatest_x_bin_global) {
          greatest_x_bin_global = greatest_x_bin_local;
        }
        if (greatest_y_bin_local > greatest_y_bin_global) {
          greatest_y_bin_global = greatest_y_bin_local;
        }
        if (greatest_z_bin_local > greatest_z_bin_global) {
          greatest_z_bin_global = greatest_z_bin_local;
        }
        histograms_[i][j][k][0].SetLineColor(tags_.colors_[k]);
        histograms_[i][j][k][1].SetLineColor(tags_.colors_[k]);
        histograms_[i][j][k][2].SetLineColor(tags_.colors_[k]);
      }
      for (int k = 0; k < tags_.components_.size(); ++k) {
        canvas.cd(j + 1 + 0);
        TString options = ((k == 0) ? "e" : "e same");
        histograms_[i][j][k][0].SetMaximum(1.1 * greatest_x_bin_global);
        histograms_[i][j][k][0].Draw(options);
        canvas.cd(j + 1 + 3);
        histograms_[i][j][k][1].SetMaximum(1.1 * greatest_y_bin_global);
        histograms_[i][j][k][1].Draw(options);
        canvas.cd(j + 1 + 6);
        histograms_[i][j][k][2].SetMaximum(1.1 * greatest_z_bin_global);
        histograms_[i][j][k][2].Draw(options);
      }
    }
    canvas.Print(tags_.species_[i] + ".eps");
  }
}

void Fit3D::saveHistograms(const TString& filename)
{
  TFile histogram_file(filename, "RECREATE");
  for (int i = 0; i < tags_.species_.size(); ++i) {
    for (int j = 0; j < tags_.signs_.size(); ++j) {
      for (int k = 0; k < tags_.components_.size(); ++k) {
        for (int l = 0; l < 3; ++l) {
          histograms_[i][j][k][l].Write();
        }
      }
    }
  }
  histogram_file.Close();
}

void Fit3D::generateModels()
{
  std::cout << "Generating models..." << endl;
  TCut pp_events("event_sign == event_sign::pp");
  TCut pn_events("event_sign == event_sign::pn");
  TCut nn_events("event_sign == event_sign::nn");
  TCut ss_events = nn_events || pp_events;
  TCut bs_events("component == component::bs");
  TCut bd_events("component == component::bd");
  TCut cc_events = bs_events || bd_events;
  TCut cw_events("component == component::cw");
  TCut ww_events("component == component::ww");
  TCut cn_events("component == component::cn");

  // Reduce data to components.
  RooDataSet* cc_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cc_events && pp_events));
  RooDataSet* cc_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cc_events && nn_events));
  RooDataSet* cc_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cc_events && pn_events));
  RooDataSet* cw_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cw_events && pp_events));
  RooDataSet* cw_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cw_events && nn_events));
  RooDataSet* cw_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cw_events && pn_events));
  RooDataSet* ww_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(ww_events && pp_events));
  RooDataSet* ww_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(ww_events && nn_events));
  RooDataSet* ww_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(ww_events && pn_events));
  RooDataSet* cn_pp_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cn_events && pp_events));
  RooDataSet* cn_nn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cn_events && nn_events));
  RooDataSet* cn_pn_data = (RooDataSet*) data_set_->reduce(
      RooFit::Cut(cn_events && pn_events));

  RooArgList observables(*x_variable_, *y_variable_, *z_variable_);
  RooDataHist cc_pp_binned_data("cc_pp_binned_data", "cc_pp_binned_data",
      observables, *cc_pp_data);
  RooDataHist cw_pp_binned_data("cw_pp_binned_data", "cw_pp_binned_data",
      observables, *cw_pp_data);
  RooDataHist ww_pp_binned_data("ww_pp_binned_data", "ww_pp_binned_data",
      observables, *ww_pp_data);
  RooDataHist cn_pp_binned_data("cn_pp_binned_data", "cn_pp_binned_data",
      observables, *cn_pp_data);
  RooDataHist cc_nn_binned_data("cc_nn_binned_data", "cc_nn_binned_data",
      observables, *cc_nn_data);
  RooDataHist cw_nn_binned_data("cw_nn_binned_data", "cw_nn_binned_data",
      observables, *cw_nn_data);
  RooDataHist ww_nn_binned_data("ww_nn_binned_data", "ww_nn_binned_data",
      observables, *ww_nn_data);
  RooDataHist cn_nn_binned_data("cn_nn_binned_data", "cn_nn_binned_data",
      observables, *cn_nn_data);

  int interpolation_order = 2;
  RooHistPdf cc_pp_pdf("cc_pp_pdf", "cc_pp_pdf",
      observables, cc_pp_binned_data, interpolation_order);
  RooHistPdf cw_pp_pdf("cw_pp_pdf", "cw_pp_pdf",
      observables, cw_pp_binned_data, interpolation_order);
  RooHistPdf ww_pp_pdf("ww_pp_pdf", "ww_pp_pdf",
      observables, ww_pp_binned_data, interpolation_order);
  RooHistPdf cn_pp_pdf("cn_pp_pdf", "cn_pp_pdf",
      observables, cn_pp_binned_data, interpolation_order);
  RooHistPdf cc_nn_pdf("cc_nn_pdf", "cc_nn_pdf",
      observables, cc_nn_binned_data, interpolation_order);
  RooHistPdf cw_nn_pdf("cw_nn_pdf", "cw_nn_pdf",
      observables, cw_nn_binned_data, interpolation_order);
  RooHistPdf ww_nn_pdf("ww_nn_pdf", "ww_nn_pdf",
      observables, ww_nn_binned_data, interpolation_order);
  RooHistPdf cn_nn_pdf("cn_nn_pdf", "cn_nn_pdf",
      observables, cn_nn_binned_data, interpolation_order);
  
  std::cout << "Finished!" << endl;
}