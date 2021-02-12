load('updated_dataset_pvc');
data(1:40,1:650000) = normal;
data(41:76,1:650000) = pvc;
Pvc = pvc(1,1:5000);
nPvc = normal(30,1:5000);
plot(Pvc,'red');
hold on
 plot(nPvc,'blue');
 xlabel('Number of Samples')
ylabel('Amplitude')
title('Normal & PVC ECG Signal');
legend('PVC Signal','Non-PVC Signal');
for i =1:numel(data(:,1))
    ff(i,:)= hist(data(i,:),500);
end
labels(1:40,1) ={'Normal'}; 
labels(41:76,1) ={'PVC'};
new = categorical(labels);

cvp = cvpartition(labels,'Holdout',0.3);
dataTrain = ff(cvp.training,:);
YTrain = labels(cvp.training,:);
dataTest = ff(cvp.test,:);
YTest = labels(cvp.test,:);
for i=1:numel(dataTrain(:,1))
B1{i,1}=dataTrain(i,:);
end
C1 = B1(:,1);
for i=1:numel(dataTest(:,1))
B2{i,1}=dataTest(i,:);
end
C2 = B2(:,1);
inputSize =1;
numHiddenUnits = 150;
numClasses = 2;
layers =[ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
    maxEpochs = 500;
    
    miniBatchSize = 256;
    options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',0, ...
    'Plots','training-progress');
trnet= trainNetwork(C1,categorical(YTrain),layers,options);