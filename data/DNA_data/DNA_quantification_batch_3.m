
%Batch 1 DNA quantification data.

n_replicates_batch_1=1;

batch_number_1=num2str(1);
    
data_batch_1=PB220307adsamplesB1A1A2modifiedforreport;

[IDs_1,Concs_1]=DNA_quant(data_batch_1,n_replicates_batch_1,batch_number_1);

hold on
title(strcat('DNA concentration for samples in experiment PB1'))
bar(IDs_1,Concs_1);
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');

%Batch 2 DNA quantification data.

n_replicates_batch_2=2;

batch_number_2=num2str(2);
    
data_batch_2=PB220307bdsamplesP1M1M2equalreplicates;

[IDs_2,Concs_2]=DNA_quant(data_batch_2,n_replicates_batch_2,batch_number_2);

hold on
title(strcat('DNA concentration for samples in experiment PB2'))
bar(IDs_2,Concs_2);
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');

%Batch 3 DNA quantification data

n_replicates_batch_3=3;

batch_number_3=num2str(3);
    
data_batch_3=PB220307cdsamplesbatch3;

[IDs_3,Concs_3]=DNA_quant(data_batch_3,n_replicates_batch_3,batch_number_3);

hold on
title(strcat('DNA concentration for samples in experiment PB3'))
bar(IDs_3,Concs_3);
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');



%Batch 4 DNA quantification data


n_replicates_batch_4=3;

batch_number_4=num2str(4);
    
data_batch_4=PB220419ddbatch4;

[IDs_4,Concs_4]=DNA_quant(data_batch_4,n_replicates_batch_4,batch_number_4);

hold on
title(strcat('DNA concentration for samples in experiment PB4'))
bar(IDs_4,Concs_4);
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');

%Batch 5 DNA quantification data

n_replicates_batch_5=3;

batch_number_5=num2str(4);
    
data_batch_5=PB220419ddbatch5;

[IDs_5,Concs_5]=DNA_quant(data_batch_5,n_replicates_batch_5,batch_number_5);

hold on
title(strcat('DNA concentration for adult samples batch 3 frozen, manually cut, and non-centrifuged with different SDS conditions'))
bar(IDs_5,Concs_5)
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');



%% Piglets only plot

n_replicates_b=3;

batch_number_b=num2str(6);
    
data_batch_b=Pigletdatabatch2and3;

[IDs_b,Concs_b]=DNA_quant(data_batch_b,n_replicates_b,batch_number_b);

piglet_table=table(IDs_b,Concs_b);

%Sorting by weight.
piglet_table_sorted=sortrows(piglet_table,2);

%Plotting the sorted table:
hold on
title(strcat('DNA concentration for all piglet samples with different conditions.'))
bar(piglet_table_sorted{1:end,1},piglet_table_sorted{1:end,2})
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');



%% Adults only plot

n_replicates_a=3;

batch_number_a=num2str(6);
    
data_batch_a=Adultdatabatch2345;

[IDs_a,Concs_a]=DNA_quant(data_batch_a,n_replicates_a,batch_number_a);

adult_table=table(IDs_a,Concs_a);

%Sorting by weight.
adult_table_sorted=sortrows(adult_table,2);

%Plotting the sorted table:
hold on
title(strcat('DNA concentration for all adult samples with different conditions. non-centrifuged is manually cut'))
bar(adult_table_sorted{1:end,1},adult_table_sorted{1:end,2})
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');

IDbrah=table2array(adult_table_sorted(:,1));
DNAbrah=table2array(adult_table_sorted(:,2));


hold on
title(strcat('DNA concentration for all adult samples with different conditions. non-centrifuged is manually cut'))
bar(IDbrah,DNAbrah)
xlabel('Samples')
ylabel('DNA concentration in ng/mg of dry weight tissue');








