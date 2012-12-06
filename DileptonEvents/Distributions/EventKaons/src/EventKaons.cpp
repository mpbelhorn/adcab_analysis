#include "EventKaons.h"

EventKaons::EventKaons(
    TString input_ntuple_file,
    TString analysis_name,
    TString x_axis_label,
    TString x_axis_unit,
    double min_x_bin_edge,
    double max_x_bin_edge)
  : DileptonEvents(input_ntuple_file, analysis_name)
{
  std::cout << "Initiliazing Fit3D class." << std::endl;
  max_event_kaons_ = 5;
  output_path_
      = "/home/matt/research/belle/adcab/analysis/DileptonEvents/Distributions"
      "/EventKaons/output/" + analysis_name + "/";
  
  boost::filesystem::path output_dir(output_path_);
  
  if (boost::filesystem::exists(output_dir)) {
    cout << "Output path exists. Overwriting directory contents!\n";
  } else {
    cout << "Output path does not exist. Creating new directory.\n";
    boost::filesystem::create_directories(output_dir);
  }
  
  // Construct the histograms.
  TString name;
  name.Clear();
  histograms_.clear();
  
  for (int species = 0; species < tags_.species_.size(); ++species) {
    EventSignHistograms event_sign_histograms;
    event_sign_histograms.clear();
    for (int dilepton_sign = 0; dilepton_sign < tags_.signs_.size();
        ++dilepton_sign) {
      EventKaonHistograms event_kaon_histograms;
      event_kaon_histograms.clear();
      for (int event_kaons = 0; event_kaons <= max_event_kaons_;
          ++event_kaons) {
        KaonSignHistograms kaon_sign_histograms;
        kaon_sign_histograms.clear();
        for (int number_of_k_minus = 0; number_of_k_minus < event_kaons + 1;
            ++number_of_k_minus) {
          ComponentHistograms component_histograms;
          component_histograms.clear();
          for (int component = 0; component < tags_.components_.size();
              ++component) {
            name.Append(tags_.species_[species]);
            name.Append("_");
            name.Append(tags_.signs_[dilepton_sign]);
            name.Append("_k");
            name += event_kaons;
            name.Append("-");
            name += number_of_k_minus;
            name.Append("_");
            name.Append(tags_.components_[component]);
            TH1D component_histogram(name, name, 
                30, min_x_bin_edge, max_x_bin_edge);
            component_histograms.push_back(component_histogram);
            name.Clear();
          }
          kaon_sign_histograms.push_back(component_histograms);
        }
        event_kaon_histograms.push_back(kaon_sign_histograms);
      }
      event_sign_histograms.push_back(event_kaon_histograms);
    }
    histograms_.push_back(event_sign_histograms);
  }
}


EventKaons::~EventKaons()
{
  // std::cout << "Starting destructor." << std::endl;
  // TODO member pointers need to be deleted, but these lines give segfault.
  // delete x_variable_;
  // delete y_variable_;
  // delete y_variable_;
}

void EventKaons::ntupleLoopCore(const int& entry_id)
{
  x_value_ = l0_pcm + l1_pcm;
  AnalysisSelectors cuts(*this);
  fillHistograms(n_kaons, n_k_min, cuts.signal_type());
}

void EventKaons::fillHistograms(
    const float& event_kaons,
    const int& number_of_k_minus,
    const int& component)
{
  // Catch bad entries. Real data will have typ_tru = 0.
  // Use typ_asn for real data.
  if (event_kaons <= max_event_kaons_) {
    histograms_[typ_asn][evt_sign + 1]
        [event_kaons][number_of_k_minus][component].Fill(x_value_);
  }
}

void EventKaons::drawHistograms()
{
  for (int species = 0; species < tags_.species_.size(); ++species) {
    for (int dilepton_sign = 0; dilepton_sign < tags_.signs_.size();
        ++dilepton_sign) {
      TString canvas_name =
          tags_.species_[species] + "_" + tags_.signs_[dilepton_sign];
      TCanvas canvas("canvas", canvas_name, max_event_kaons_ * 500,
          max_event_kaons_ * 500);
      canvas.Divide(max_event_kaons_, max_event_kaons_);
      double canvas_maximum_scale = 0;
      for (int event_kaons = 0; event_kaons <= max_event_kaons_;
          ++event_kaons) {
        for (int number_of_k_minus = 0; number_of_k_minus < event_kaons + 1;
            ++number_of_k_minus) {
          TH1D* total = (TH1D*) histograms_.
              at(species)[dilepton_sign][event_kaons][number_of_k_minus][0].
              Clone("Total");
          for (int component = 1; component < tags_.components_.size();
              ++component) {
            total->Add(&histograms_[species][dilepton_sign][event_kaons]
                [number_of_k_minus][component]);
          }
          if (event_kaons == 0) {
            canvas_maximum_scale = total->GetMaximum();
          } else {
            total->SetMaximum(canvas_maximum_scale);
          }
          canvas.cd(number_of_k_minus + 1 + event_kaons * max_event_kaons_);
          total->SetLineColor(1);
          total->SetMinimum(0);
          total->Draw("e");
            
          for (int component = 0; component < tags_.components_.size();
              ++component) {
            TString options = "e same";
            histograms_[species][dilepton_sign][event_kaons][number_of_k_minus]
                [component].SetLineColor(tags_.colors_[component]);
            histograms_[species][dilepton_sign][event_kaons][number_of_k_minus]
                [component].Draw(options);
          }
        }
      }
      canvas.Print(output_path_ + canvas_name + ".eps");
    }
  }
}


void EventKaons::saveHistograms(const TString& filename)
{
  TFile histogram_file(filename, "RECREATE");
  for (int species = 0; species < tags_.species_.size(); ++species) {
    for (int dilepton_sign = 0; dilepton_sign < tags_.signs_.size();
        ++dilepton_sign) {
      for (int event_kaons = 0; event_kaons <= max_event_kaons_;
          ++event_kaons) {
        for (int number_of_k_minus = 0; number_of_k_minus < event_kaons + 1;
            ++number_of_k_minus) {
          for (int component = 1; component < tags_.components_.size();
              ++component) {
            histograms_[species][dilepton_sign][event_kaons][number_of_k_minus]
                [component].Write();
          }
        }
      }
    }
  }
  histogram_file.Close();
}
