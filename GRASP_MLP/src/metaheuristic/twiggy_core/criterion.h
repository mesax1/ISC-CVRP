#ifndef TWIGY_CORE_CRITERION_H_
#define TWIGY_CORE_CRITERION_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "utils.h"

enum ImpurityMeasure { gini, entropy };

class Criterion {
 public:
  const std::vector<std::vector<float>> *feature_data_;
  const std::vector<int> *label_data_;
  ImpurityMeasure impurity_measure_;
  const std::vector<SampleData> *sample_map_ptr_;
  float (*impurity_fn_)(std::vector<int> &, int &, int &);
  int start_;
  int end_;
  int pos_;
  int n_labels_;
  std::vector<int> label_freqs_total_;
  std::vector<int> label_freqs_left_;
  std::vector<int> label_freqs_right_;
  int n_samples_;
  int n_samples_left_;
  int n_samples_right_;

  Criterion(ImpurityMeasure impurity_measure, int n_labels,
            const std::vector<std::vector<float>> *feature_data,
            const std::vector<int> *label_data);
  Criterion() {}
  void SetData(const std::vector<SampleData> *sample_map_ptr);
  void SetNodeLimits(int start, int end);
  void UpdateSplitPos(int new_pos);
  float ImpurityImprovement();
  float NodeImpurity();
  void ChildrenImpurities(float &impurity_left, float &impurity_right);
  void ResetStats();
  static float GiniCoefficient(std::vector<int> &label_freqs, int &n_samples,
                                int &n_labels);
  static float Entropy(std::vector<int> &label_freqs, int &n_samples,
                        int &n_labels);
};

#endif  // TWIGY_CORE_CRITERION_H_
