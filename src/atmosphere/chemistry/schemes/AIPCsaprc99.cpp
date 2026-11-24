//************************************************************************
// This computer software was prepared by University of Science and Technology of China (USTC)
// NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
// EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
//
// Module to compute AIPC_saprc99 subroutine for chem(:,:,:,:).
// * Author: Zihan Xia, ChenJin
// * Email: xiazh@mail.ustc.edu.cn
// Last update: zhxia 09/2024
//
// Contact:
// Chun Zhao, PhD
// Professor
// University of Science and Technology of China
// Hefei, Anhui, 230026
// Email: ChunZhao@ustc.edu.cn
//
//
// Additional information:
// aemol.ustc.edu.cn: resources
//
//************************************************************************

#include <torch/script.h> // One-stop header.
#include "torch/torch.h"

#include <iostream>
#include <memory>

extern "C" void AIPC_saprc99(double *c_cin, double *c_cout, int batch_size);

void AIPC_saprc99(double *c_cin, double *c_cout, int batch_size) 
{
  const char * model_path;
  model_path = "AIPC_saprc99.pt";

  torch::jit::script::Module module;
  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    module = torch::jit::load(model_path);
  }
  catch (const c10::Error& e) {
    std::cerr << "error loading the model\n";
    return;
  }
 
  auto cinput_ = torch::from_blob(c_cin, {batch_size, 91}, torch::kDouble).unsqueeze(1); //(batch_size,channel=1,91 features)
  cinput_ = cinput_.to(torch::kFloat);
 
  // Create a vector of inputs.
  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(cinput_);   

// Execute the model and turn its output into a tensor.
  at::Tensor output = module.forward(inputs).toTensor();
  output = output.to(torch::kDouble);

// Copy the data from the output tensor to c_cout.
  std::memcpy(c_cout, output.data_ptr<double>(), 74 * batch_size * sizeof(double));
}
