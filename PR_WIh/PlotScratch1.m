data=load('PRwIhTTFSAMPA4curveslikefig17a.mat');
figure()
for i=1:size(data.ghVdsOutTTFS,1)
    for j=1:size(data.ghVdsOutTTFS,2)
        tmpVdsOut=squeeze(data.ghVdsOutTTFS(i,j,:,1));
        tmpTTFS=squeeze(data.ghVdsOutTTFS(i,j,:,2));
        plot(tmpVdsOut,tmpTTFS)
        hold on;
    end
end
        
