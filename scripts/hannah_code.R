#running Hannah's code
#have created an IVM model with x3 SEI compartments, depending on whether the mosquito gets IVM from human/cattle/none at all
require(ICDMM)
require(tidyverse)
require(cowplot)

#base model
out <- run_model(model = "odin_model_IVM_SMChet",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   time = 730)
