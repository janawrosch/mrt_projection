% mrt_projection_fertig_mrt.m
% Program to project diffusion MRI data to a two-dimensional map
% Saved data that is being read:
% load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\mrt_2sigma_karte.mat')
% 	2-sigma-reference map of what diffusion levels are expected in each pixel in a healthy brain
% 	Adjust path also in line 142
% load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\cmap_mrt.mat');
% 	Color map for two-dimensional map of the brain
% 	Adjust path also in line 194
% path = ('V:\AG_Neurophotonik\Projekte\Projektion\Daten\all_cases\');
% 	Path in which diffusion MRI data is saved as groups of .dcm-files
% Input:
% startcase: First casenumber that is being projected
% ncases: Last casenumber that is being projected
% (If only one case is projected give startcase and ncases the value of that casenumber)
% loaddata:	Set to 1 if the .dcm-files are new and need to be projected for the first time
% 		Set to 0 if the projected data is already stored in a fertiges_img2.mat-file
% Index-Table: If a stack of .dcm-files for a certain case is not numbered 1 to 24 adjust the MATLAB code here to fit the available data.
% Output:
% stroke_brain_ratio_MRI: Number of stroke voxels divided by the number of all brain voxels in the diffusion MRI data
% stroke_brain_ratio_map: Number of stroke pixels divided by the number of all brain pixels in the two-dimensional stroke map
% delta_ratio: absolute difference between stroke_brain_ratio_MRI and stroke_brain_ratio_map. This is a measure of projection quality. If it is too high (higher than delta_ratio_threshold in territory_recognition.m) this projection failed and will be excluded from further analysis.
% fertiges_img2.mat: projected data will be saved here and can be read in directly next time by setting loaddata to 0.






clear all
close all
set(0, 'DefaulttextInterpreter', 'none')

 %%%%%%% Abgespeicherte Daten auf die zugegriffen wird    %%%%%%%%%%%%%%%
 load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\mrt_1sigma_karte.mat') %2sigma Mittelwert-Karte; steht nochmal in Zeile 142
 load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\cmap_mrt.mat');     %dazu passende Colormap; steht nochmal in Zeile 194
 path = ('V:\AG_Neurophotonik\Projekte\Projektion\Daten\all_cases\'); %Pfad in dem die .tif Bilder zu finden sind und in dem die Reusltate abgespeichert werden        

 %%%%%%% EINGABEN %%%%%%%%%%%%%%%%%%%%%
startcase=1; % Patient mit dem angefangen werden soll
ncases=105; %letzter Patient der ausgewertet werden soll




h=waitbar(ncases,'Projection...');
for casenumber=startcase:ncases
  waitbar(casenumber/ncases)  
        filename = sprintf('%i_', casenumber); %Dateiname der .tif-Bilder ohne Nummerierung (soviele Nullen stehenlassen wie beim höchsten Index noch vor der Zahl stehen)
        savename=sprintf('%i_fertiges_img2.mat', casenumber);
        
        
% INDEX-TABLE        
         startstop=[1 24]; %Nummerierung der einzulesenden .tif Bilder (Index des ersten und letzten Bildes)
        if casenumber==17
            startstop=[1 22]; %Nummerierung der einzulesenden .tif Bilder (Index des ersten und letzten Bildes)      
        elseif casenumber==20
            startstop=[1 26]; %Nummerierung der einzulesenden .tif Bilder (Index des ersten und letzten Bildes)              
        elseif casenumber==47
            startstop=[1 17]; %Nummerierung der einzulesenden .tif Bilder (Index des ersten und letzten Bildes)                       
        elseif casenumber == 63
            startstop=[1 21]; %Nummerierung der einzulesenden .tif Bilder (Index des ersten und letzten Bildes)            
        elseif casenumber == 69
            startstop=[1 23]; %Nummerierung der einzulesenden .tif Bilder (Index des ersten und letzten Bildes)      
    end
    

     
%%%%%%% ENDE DER EINGABE %%%%%%%%%%%%%%%%%%     
  



loaddata=0; % 1 wenn  Daten neu eingelesen werden; 0 wenn die MRT Daten schon in img2 gespeichert sind
tic  
if loaddata==1
%% Laden der Daten
    [stack1]=limg(path,filename,'.dcm', startstop,1);
    stack1 = double(stack1);
%     data=stack1;
    data = double(flipdim(stack1, 3));
   

    
 
        
    cchoice = 'global';
    
    % Mittelpunkt, Orientierung und Threshold festlegen
    dynimghist(data, cchoice);
    
    global dynimghist_point dynimghist_point2 dynimghist_threshold   strokethreshold
%     dynimghist_point = [65.0050   73.1924   11.5000];
%     dynimghist_point2 = [64.0000   14.0000   11.5000];                  %Default für Schmidt, Sigrid
%     dynimghist_threshold = 89.6100;
%     dynimghist_colormapthreshold = 250.2900;
    
    center = dynimghist_point;
    angle = rad2deg(cart2pol(dynimghist_point2(1) - center(1), dynimghist_point2(2) - center(2)));
    threshold = dynimghist_threshold;

    %% Kugelkoordinaten berechnen
    layer = 2;
    ndatapoints = sum(sum(sum(data > threshold))); %ndatapoints=numel(find(data>threshold));
    nstroke_MRI=sum(sum(sum(data>strokethreshold))); % number of stroke pixels in MRI
    datapoints = zeros(ndatapoints, 3);
    
    count = 0;
    for y = 1:size(data, 1)
        for x = 1:size(data, 2)
            for z = 1:size(data, 3) %für jeden MRT Pixel
                if (data(y, x, z) > threshold) %wenn der Pixel im MRT> threshold
                    count = count + 1;
                    [datapoints(count, 2), datapoints(count, 1)] = cart2sph(x - center(1), y - center(2), (z - center(3))*layer); %Zentrum des Koordinatensystems ist im festgelegten Punkt. Cart-Koordinaten (vom Zentrum aus) werden 
                    % umgerechnet in Kugelkoordinaten, dabei interessieren uns nur theta und phi, da r für alle Punkt gleich ist und sie so auf die Kugeloberfläche projiziert werden
                    a=datapoints;
                    datapoints(count, 1) = mod(rad2deg(datapoints(count, 1)) + 90, 180) - 90; %Umrechnung der Kugel(oberflächen)koordinaten von Bogemaß in Gradmaß
                    datapoints(count, 2) = mod(rad2deg(datapoints(count, 2)) - angle, 360);
                    datapoints(count, 3) = data(y, x, z);
                end
            end
        end
    end
    
    %% Projektion der Kugel in die Ebene
    mstruct = defaultm('mollweid');
    mstruct.origin = [0, 180, 0];
    mstruct.maplatlimit = [-90, 90];
    mstruct.maplonlimit = [0, 360];
    mstruct = defaultm(mstruct);
    
    [x, y] = mfwdtran(mstruct, datapoints(:, 1), datapoints(:, 2));
    projdatapoints = [x, y, datapoints(:, 3)];

    %% Bild erzeugen
    x = x - min(x); x = 360*x/max(x);
    y = y - min(y); y = 180*y/max(y);
    img = zeros(181, 361);
    
    for k = 1:size(projdatapoints, 1)
        if (img(round(y(k)) + 1, round(x(k)) + 1) < projdatapoints(k, 3))
            img(round(y(k)) + 1, round(x(k)) + 1) = projdatapoints(k, 3);
        end
    end
    
    img = fliplr(img);
    
    %% Interpolation des Bildes
    mwell = load('mollweidellipsoid');
    mwell = mwell.mollweidellipsoid;
    
    img2 = zeros(181, 361);
    imgt = zeros(185, 365);
    imgt(3:183, 3:363) = img;
    for i = 3:183
        for j = 3:363
            img2(i-2, j-2) = mean(mean(imgt((i-2):(i+2), (j-2):(j+2))));
        end
    end
    
    img2 = img2.*mwell;
    
    
save([path, savename], 'img2', 'dynimghist_colormapthreshold');
   
    elseif loaddata== 0
        load([path, savename]);
    end
    
%% Abzug der 2sigma-Karte und Bild
load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\mrt_1sigma_karte.mat') %2sigma Mittelwert-Karte
%Normierung der Daten
  
imgrange = [1, 99]; imgrange(3) = imgrange(2) - imgrange(1) + 1;


    img2max=max(max(img2(:,:)));
    img2min=min(min(img2(:,:)));
    img2(:,:)=imgrange(3)*(img2(:,:)-img2min)/(img2max - img2min) + imgrange(1);

%Abzug der 2sigma_karte

fertiges_img2=zeros(size(img2));

nbrain_map=0;
nstroke_map=0;
for i=1:size(img2,1)
    for j=1:size(img2,2)
        if img2(i,j)==1
            fertiges_img2(i,j)=0;
        elseif img2(i,j)-mrt_2sigma_karte(i,j) <=0
            nbrain_map=nbrain_map+1;
            fertiges_img2(i,j)=1;
        elseif img2(i,j)-mrt_2sigma_karte(i,j) > 0
            fertiges_img2(i,j)= 50+ img2(i,j)-mrt_2sigma_karte(i,j);
            nstroke_map=nstroke_map+1;
        end
    end
end

%Normierung für Colormap
imgrange = [1, 99]; imgrange(3) = imgrange(2) - imgrange(1) + 1;

    img2max=max(fertiges_img2(fertiges_img2>1));
    img2min=min(fertiges_img2(fertiges_img2>1));

    for i=1:size(fertiges_img2,1)
        for j=1:size(fertiges_img2,2)    
            if fertiges_img2(i,j)>1
    fertiges_img2(i,j)=imgrange(3)*(fertiges_img2(i,j)-img2min)/(img2max - img2min) + imgrange(1);
            end
        end
    end
    
    stroke_brain_ratio_MRI=nstroke_MRI/ndatapoints;
    stroke_brain_ratio_map=nstroke_map/nbrain_map;
    delta_ratio=abs(stroke_brain_ratio_map-stroke_brain_ratio_MRI);
    
keep_ratios(casenumber,1)=delta_ratio;

    %% Bilder anzeigen
    % MRT projection
 load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\cmap_mrt.mat');     %dazu passende Colormap    
    fertiges_img2=flipdim(fertiges_img2,2);
    figure('name', 'MRT projection');
    imagesc(-180:180, -90:90, fertiges_img2); 
    axis xy;
    colormap(cmap);
    title(sprintf('MRT-Projektion: %i \n Stroke/Brain Ratio MRI: %i\n Stroke/Braion Ratio map: %i \n Difference in Ratio: %i', casenumber, stroke_brain_ratio_MRI, stroke_brain_ratio_map, delta_ratio));

  time=toc;
 path = ('V:\AG_Neurophotonik\Projekte\Projektion\Daten\all_cases\'); %Pfad in dem die .tif Bilder zu finden sind und in dem die Reusltate abgespeichert werden        
 save([path, savename]);

end

close(findall(0,'Tag','TMWWaitbar'));
