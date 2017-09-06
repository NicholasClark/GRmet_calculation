files = dir('*csv');
files = {files.name};
t_all = table;
for i=1:3
    t_ = readtable(files{i});
    t_ = [t_(:,[2 4 5  10 11 12]) ...
        table(i*ones(height(t_),1), 'variablenames', {'replicate'})];
    t_.Properties.VariableNames = {
        'cell_line'
        'small_molecule'
        'concentration'
        'cell_count__time0'
        'cell_count'
        'cell_count__ctrl'
        'replicate'};
    t_all = [t_all; t_];
end

writetable(t_all, 'LJP_data_replicates.tsv', ...
        'filetype','text' , 'delimiter', '\t')
%table2tsv(t_all, 'LJP_data_replicates.tsv')

%%

addpath('../MATLAB')


[t_GRvalues, t_GRmetrics] = GRmetrics('', 'LJP_data_replicates.tsv');

writetable(t_GRmetrics, 'LJP_GRmetrics_replicates_NC2.tsv', ...
        'filetype','text' , 'delimiter', '\t')
writetable(t_GRvalues, 'LJP_GRvalues_replicates_NC2.tsv', ...
        'filetype','text' , 'delimiter', '\t')
    
%table2tsv(t_GRvalues, 'LJP_GRvalues_replicates_NC.tsv')
%table2tsv(t_GRmetrics, 'LJP_GRmetrics_replicates_NC.tsv')

[t_GRvalues, t_GRmetrics] = GRmetrics('', 'LJP_data_replicates.tsv', ...
    'collapseKeys', 'replicate');
%table2tsv(t_GRvalues, 'LJP_GRvalues_merged_NC.tsv')
%table2tsv(t_GRmetrics, 'LJP_GRmetrics_merged_NC.tsv')
writetable(t_GRmetrics, 'LJP_GRmetrics_merged_NC2.tsv', ...
        'filetype','text' , 'delimiter', '\t')
writetable(t_GRvalues, 'LJP_GRvalues_merged_NC2.tsv', ...
        'filetype','text' , 'delimiter', '\t')