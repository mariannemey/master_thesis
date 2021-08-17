clear

%Path to images and extracting image file infromation 
data_path = 'C:\Users\maria_7lofnoj\OneDrive\Desktop\New folder\NENT181';
addpath(data_path)
file_information = fullfile(data_path, '*.jpg');
jpegFiles = dir(file_information);

%Initiating results table
rows = size(jpegFiles,1);
table_size = [rows 2];
table_variable_types = {'string','double'};
table_variable_names = {'Image', 'Percentage'};
results_table = table('Size', table_size,'VariableTypes', table_variable_types, 'VariableNames', table_variable_names);

tic
%Analysing images and filling results table
for k = 1:length(jpegFiles)
    %Loading image
    baseFileName = jpegFiles(k).name;
    fullFileName = fullfile(data_path, baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
    image =  imread(fullFileName);
    image_1_meta = imfinfo(fullFileName); %MM - reads meta data - data about the data
    %imtool(image,[])
    
    %Unstacking image
    ch1_image = image(:,:,1);
    ch2_image = image(:,:,2);
    ch3_image = image(:,:,3);
    %imtool(ch1_image,[])
    %imtool(ch2_image,[])
    %imtool(ch3_image,[])
    
    %Creating cell mask 
    cell_mask = ch2_image > 55; % Threshold
    %imtool(cell_mask, [])
    cell_mask = ~cell_mask;
    %imtool(cell_mask, [])
    
    %Analysing mask 
    area = sum(cell_mask(:));
    total_area_possible = size(image, 1) * size(image, 2);
    percentage = (area/total_area_possible)*100;
    
    %Filling results table 
    results_table.Image(k) = string(baseFileName);
    results_table.Percentage(k) = percentage; 
end
toc

%Exporting results table 
addpath(data_path)
%.txt
writetable(results_table,'results_table.txt');
%.xls
writetable(results_table,'results_MEMT181.xls');
