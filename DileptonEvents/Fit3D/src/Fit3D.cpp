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
    TString species_string = event_species->GetName();
    HistogramListOverEventSigns event_sign_histograms;
    while((event_sign = (RooCatType*) event_sign_iterator->Next())) {
      TString sign_string = event_sign->GetName();
      HistogramListOverComponents component_histograms;
      component_histograms.clear();
      while((component = (RooCatType*) component_iterator->Next())) {
        HistogramListOverVariables variable_histograms;
        variable_histograms.clear();
        TString component_string = component->GetName();
        name.Append(species_string);
        name.Append("_");
        name.Append(sign_string);
        name.Append("_");
        name.Append(component_string);
        
        TaggedHistogram x;
        TaggedHistogram y;
        TaggedHistogram z;
        
        x.event_species = species_string;
        x.event_sign = sign_string;
        x.component = component_string;
        x.variable = x_axis_label;
        x.histogram = TH1F(name + "_x", name + "_x", 100, min_x_bin_edge, max_x_bin_edge);
        variable_histograms.push_back(x);
        y.event_species = species_string;
        y.event_sign = sign_string;
        y.component = component_string;
        y.variable = y_axis_label;
        y.histogram = TH1F(name + "_y", name + "_y", 100, min_y_bin_edge, max_y_bin_edge);
        variable_histograms.push_back(y);
        z.event_species = species_string;
        z.event_sign = sign_string;
        z.component = component_string;
        z.variable = z_axis_label;
        z.histogram = TH1F(name + "_z", name + "_z", 100, min_z_bin_edge, max_z_bin_edge);
        variable_histograms.push_back(z);
        
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
//  if (typ_tru == 0 || component == 0) {
    // std::cout << "Bad entry found with species " << typ_tru
    //           << " and component " << component << endl;
 // } else {
      histograms_[typ_tru][evt_sign + 1][component][0].histogram.Fill(x_value_);
      histograms_[typ_tru][evt_sign + 1][component][1].histogram.Fill(y_value_);
      histograms_[typ_tru][evt_sign + 1][component][2].histogram.Fill(z_value_);
  //}
}

void Fit3D::drawHistograms()
{
  HistogramListOverVariables::iterator variable_histograms;
  HistogramListOverComponents::iterator component_histograms;
  HistogramListOverEventSigns::iterator sign_histograms;
  HistogramListOverEventSpecies::iterator histograms;
  for (histograms = histograms_.begin();
      histograms < histograms_.end(); histograms++) {
    TString species_string = (*histograms).front().
        front().front().event_species;
    TCanvas canvas("canvas", species_string, 1200, 1200);
    canvas.Divide(3, 3);
    int canvas_column(0);
    for (sign_histograms = (*histograms).begin();
        sign_histograms < (*histograms).end(); sign_histograms++) {
      canvas_column++;
      int component_number = 0;
      double greatest_x_bin_global(0);
      double greatest_y_bin_global(0);
      double greatest_z_bin_global(0);
      for (component_histograms = (*sign_histograms).begin();
          component_histograms < (*sign_histograms).end();
          component_histograms++) {
        TString component_string = (*component_histograms)[0].
            component;
        double greatest_x_bin_local = (*component_histograms)[0].
            histogram.GetMaximum();
        double greatest_y_bin_local = (*component_histograms)[1].
            histogram.GetMaximum();
        double greatest_z_bin_local = (*component_histograms)[2].
            histogram.GetMaximum();
        if (greatest_x_bin_local > greatest_x_bin_global) {
          greatest_x_bin_global = greatest_x_bin_local;
        }
        if (greatest_y_bin_local > greatest_y_bin_global) {
          greatest_y_bin_global = greatest_y_bin_local;
        }
        if (greatest_z_bin_local > greatest_z_bin_global) {
          greatest_z_bin_global = greatest_z_bin_local;
        }
      }
      for (component_histograms = (*sign_histograms).begin();
          component_histograms < (*sign_histograms).end();
          component_histograms++) {
        component_number++;
        canvas.cd(canvas_column + 0);
        TString options = ((component_number == 1) ? "e" : "e same");
        (*component_histograms)[0].histogram.
            SetMaximum(1.1 * greatest_x_bin_global);
        (*component_histograms)[0].histogram.Draw(options);
        canvas.cd(canvas_column + 3);
        (*component_histograms)[1].histogram.
            SetMaximum(1.1 * greatest_y_bin_global);
        (*component_histograms)[1].histogram.Draw(options);
        canvas.cd(canvas_column + 6);
        (*component_histograms)[2].histogram.
            SetMaximum(1.1 * greatest_z_bin_global);
        (*component_histograms)[2].histogram.Draw(options);
      }
    }
    canvas.Print(species_string + ".eps");
  }
}