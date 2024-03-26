clear;
load('inputdata.mat')

num_features = 73;
num_classes = 3;
output_dimensions = 10;
hiddenLayerSize = 300;

layers = [
    featureInputLayer(num_features) % Feature input layer for 1D data
    fullyConnectedLayer(hiddenLayerSize, WeightsInitializer="he")
    eluLayer % ELU activation
    fullyConnectedLayer(hiddenLayerSize, WeightsInitializer="he")
    eluLayer % ELU activation
    fullyConnectedLayer(hiddenLayerSize, WeightsInitializer="he")
    eluLayer % ELU activation
    fullyConnectedLayer(hiddenLayerSize, WeightsInitializer="he")
    eluLayer % ELU activation
    fullyConnectedLayer(hiddenLayerSize, WeightsInitializer="he")
    eluLayer % ELU activation
    fullyConnectedLayer(num_classes * output_dimensions)
    sigmoidLayer % Sigmoid activation for multi-label classification
    regressionLayer
];
%lgraph = layerGraph(layers);
%plot(lgraph);


% Generate random input data and label matrices
num_samples = 12000;
input_data = features;
output_labels = labels;
idx = randperm(12000,2000);
XValidation = input_data(idx,:);
input_data(idx,:) = [];
YValidation = output_labels(idx,:);
output_labels(idx,:) = [];

% Define options for training with progress plot
options = trainingOptions('adam', ...
    'MaxEpochs', 25, ...
    'MiniBatchSize', 64, ...
    'InitialLearnRate', 0.001, ...
    'ValidationData',{XValidation,YValidation}, ...
    'Plots', 'training-progress'); % Enable progress plot



% Create and train the neural network
dnn = trainNetwork(input_data, output_labels, layers, options);

new_data = testFeatures;
new_predictions = predict(dnn, new_data);

outputdata_file = matfile('outputdata.mat','Writable',true);
outputdata_file.pred = new_predictions;