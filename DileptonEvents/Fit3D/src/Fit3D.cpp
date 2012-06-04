#include "Fit3D.h"
#include "EventSelectors.h"
#include <iostream>


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
  
  component_->printMultiline(std::cout, 2);
  std::cout << std::endl;
  
  // Construct the histograms.
  RooCatType *component;
  RooCatType *event_sign;
  RooCatType *event_species;
  TIterator* component_iterator = component_->typeIterator();
  TIterator* event_sign_iterator = event_sign_->typeIterator();
  TIterator* event_species_iterator = event_species_->typeIterator();
  
  event_species_iterator->Reset();
  event_sign_iterator->Reset();
  component_iterator->Reset();
  TString name;
  name.Clear();
  histograms_.clear();
  while((event_species = (RooCatType*) event_species_iterator->Next())) {
    vector<vector<vector<TH1D> > > event_sign_histograms;
    while((event_sign = (RooCatType*) event_sign_iterator->Next())) {
      vector<vector<TH1D> > component_histograms;
      component_histograms.clear();
      while((component = (RooCatType*) component_iterator->Next())) {
        vector<TH1D> variable_histograms;
        variable_histograms.clear();
        name.Append(event_species->GetName());
        name.Append("_");
        name.Append(event_sign->GetName());
        name.Append("_");
        name.Append(component->GetName());
        // std::cout << name << std::endl;
        variable_histograms.push_back(
            TH1D(name + "_x", name + "_x", 100, min_x_bin_edge, max_x_bin_edge)
        );
        variable_histograms.push_back(
            TH1D(name + "_y", name + "_y", 100, min_y_bin_edge, max_y_bin_edge)
        );
        variable_histograms.push_back(
            TH1D(name + "_z", name + "_z", 100, min_z_bin_edge, max_z_bin_edge)
        );
        component_histograms.push_back(variable_histograms);
        name.Clear();
      }
      event_sign_histograms.push_back(component_histograms);
      component_iterator->Reset();
    }
    histograms_.push_back(event_sign_histograms);
    event_sign_iterator->Reset();
  }
  event_species_iterator->Reset();
}


Fit3D::~Fit3D()
{
  std::cout << "Starting destructor." << std::endl;
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
  fillHistograms(cuts.tru.type(), cuts.event_sign(), cuts.signal_type());
  
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

void Fit3D::fillHistograms(
      const char& species,
      const char& sign,
      const int& component)
{
  // Catch bad entries.
  if (species == 0 || component == 0) {
    std::cout << "Bad entry found with species " << species
              << " and component " << component << endl;
  } else {
    double values[] = {x_value_, y_value_, z_value_};
    for (int variable = 0; variable < 3; variable++) {
      histograms_[0][sign + 1][0][variable].Fill(values[variable]);
      histograms_[species][sign + 1][0][variable].Fill(values[variable]);
      histograms_[species][sign + 1][component][variable].Fill(values[variable]);
    }
  }
}

void Fit3D::drawHistograms()
{
  // For each species: (x4)
    // For each event sign: (x3)
      // For each component: (x6+1)
        // For each variable: (x3)
          // Create a histogram.
        // For each combination of two variables: (x3)
          // Create a 2D histogram.
  TH1* hist = data_set_->createHistogram("name", *x_variable_, RooFit::AutoBinning(100, 1));
  
  TCanvas canvas("canvas", "Data Components", 1200, 1600);
  canvas.Divide(1, 2);
  canvas.cd(1);
  hist->Draw();
  canvas.Print("test.eps");
}