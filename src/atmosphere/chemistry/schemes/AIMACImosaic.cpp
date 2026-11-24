//************************************************************************
// This computer software was prepared by University of Science and Technology of China (USTC)
// NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
// EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
//
// Module to compute AIMACI_mosaic subroutine for chem(:,:,:,:).
// * Author: Zihan Xia
// * Email: xiazh@mail.ustc.edu.cn
// Last update: zhxia 10/2024
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

extern "C" void AIMACI_mosaic(double *c_cin, double *c_cout, int batch_size);

void AIMACI_mosaic(double *c_cin, double *c_cout, int batch_size) 
{
  const char * model_path;
  model_path = "AIMACI_mosaic.pt";

  torch::jit::script::Module module;
  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    module = torch::jit::load(model_path);
  }
  catch (const c10::Error& e) {
    std::cerr << "error loading the model\n";
    return;
  }
 
  auto cinput_ = torch::from_blob(c_cin, {batch_size, 67}, torch::kDouble).to(torch::kFloat); //(batch_size,67 features)
  //std::cout << "cinput_ shape: " << cinput_.sizes() << std::endl;
 
  // Create a vector of inputs.
  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(cinput_);   

// Execute the model and turn its output into a tensor.
  at::Tensor output = module.forward(inputs).toTensor();
  output = output.to(torch::kDouble); 

// Copy the data from the output tensor to c_cout.
  std::memcpy(c_cout, output.data_ptr<double>(), 61 * batch_size * sizeof(double));
}
