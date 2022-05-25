function [SampleIDcat, DNA_conc] = DNA_quant(data,n_replicates,batch_number)
%Takes a table containing the DNA data (Do not include columns past nuclueic acid concentration when importing data,
%table containing sample weights (must be in alphabetic and lowest number first order), number of replicates,
%and batch number as input and provides plots of DNA dry weight concentration
%for each sample as output. 

%Sorting data
Sorted_data=sortrows(data,2);

%Extracting columns containing SampleID, DNA concentration, and weights from
%data.
true_data=Sorted_data(:,[2,5,13]);

%Isolating the Sample ID column and only looking at one replicate.
SampleID=table2array(true_data(1:n_replicates:size(true_data,1)-n_replicates+1,1));

%Creating SampleID as a categorical for output.

SampleIDcat=categorical(SampleID);

%Isolating the weights column and only looking at one replicate.
weights=table2array(true_data(1:n_replicates:size(true_data,1)-n_replicates+1,3));

%Creating a table containing SampleID and weights.
sample_weights=table(SampleID,weights);

%Defining counting variable.
x=1;
%Calculating means of replicates of samples.
for i=1:n_replicates:size(true_data,1)-n_replicates+1
    mean_vals(x)=mean(table2array(true_data(i:i+n_replicates-1,2)));
    x=x+1;
end

%Transposing the means.
mean_vals=mean_vals.';

%Calculating DNA concentration in ng/mg of dry weight tissue.
DNA_conc=(mean_vals.*400)./table2array(sample_weights(:,2));

%Converting the batch number to a string.
batch_number=num2str(batch_number);

%Plotting the DNA concentrations of samples as a bar plot.
%hold on
%title(strcat('DNA concentration for samples in batch',{' '},batch_number))
%bar(SampleID,DNA_conc);
%xlabel('Samples')
%ylabel('DNA concentration in ng/mg of dry weight tissue');



