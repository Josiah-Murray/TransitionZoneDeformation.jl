using Images
using FileIO

pixelWidth = 100
imgSeqFolderFilePath = "Arbitrary Transition Zones/Experiments/Image Sequence/"

ArrayToConvert = copy(deformations)

#Preprocess to normalise

ArrayToConvert = ArrayToConvert .- minimum(ArrayToConvert)
ArrayToConvert = ArrayToConvert/(maximum(ArrayToConvert))


num_tValues = size(ArrayToConvert)[1]
num_tValues = 12
num_xValues = size(ArrayToConvert)[2]


for ti in 10:num_tValues
  filepath = imgSeqFolderFilePath * "deformation_t_" * string(ti) * ".png"
  img = ArrayToConvert[ti,:]
  for i in 1:pixelWidth-1
    img = [img ArrayToConvert]

  end
  save(filepath, img)
end
