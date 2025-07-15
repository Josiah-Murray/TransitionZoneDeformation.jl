using Images
using FileIO

pixelWidth = 100
imgSeqFolderFilePath = "Arbitrary Transition Zones/Experiments/Image Sequence/"

ArrayToConvert = copy(deformations)

#Preprocess to normalise

ArrayToConvert = ArrayToConvert .- minimum(ArrayToConvert)
ArrayToConvert = ArrayToConvert/(maximum(ArrayToConvert))


num_tValues = size(ArrayToConvert)[1]
num_xValues = size(ArrayToConvert)[2]

img = zeros(Float64, num_xValues, pixelWidth)

for ti in 1:num_tValues
  filepath = imgSeqFolderFilePath * "deformation_t_" * string(ti) * ".png"
  img = zeros(Float64, num_xValues, pixelWidth)
  for xi in 1:num_xValues
    for i in 1:pixelWidth-1
      img[xi,i] = ArrayToConvert[ti,xi]
    end
  end
  save(filepath, img)
end
