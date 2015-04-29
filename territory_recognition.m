

clear all 
load('V:\AG_Neurophotonik\Projekte\Projektion\Daten\all_cases\projdata.mat')
 load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\allregions.mat') 
  load('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\results.mat')
% strokepx_per_region_cut_offs=[0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 0.99];
% delta_ratio_cut_offs=[0.01 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2];
% region_covered_cut_offs=[0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 0.99];

strokepx_per_region_cut_offs=[ 1 ];
delta_ratio_cut_offs=[1 ];
region_covered_cut_offs=[  1];


numbcases=105;
count=1;
thalamus_cases=[ 26 34 17 40 49 93 52 65 66 67 78 7 21 94];
%     thalamus_cases=[3 9 12 19 22 29 32 33 54 57 58 59 60 64 71 73 76 80 26 34 17 40 49 93 52 65 66 67 78 7 21 94];
% ROC_results=zeros(size(strokepx_per_region_cut_offs,2),21);
 




for a=1:size(strokepx_per_region_cut_offs,2)
    for b=1:size(region_covered_cut_offs,2)
        for c=1:size(region_covered_cut_offs,2)
            for d=1:size(region_covered_cut_offs,2)
                for e=1:size(region_covered_cut_offs,2)
                    for f=1:size(region_covered_cut_offs,2)
                        for g=1:size(delta_ratio_cut_offs,2)
    
    
    
    
    
results(:,:,1)=zeros(numbcases,6);
for casenumber=1:numbcases
 
    %Für jeden Fall neu berechnen
results(casenumber,6,1)=projdata.delta_ratio{casenumber,1};
stroke_map=find(projdata.fertiges_img2{casenumber,1}>1); %Indizes der Strokepixel


strokepx_in_ant=sum(ismember(stroke_map(:,1),stroke_ant(:,1)));
strokepx_in_vMSG=sum(ismember(stroke_map(:,1),stroke_vMSG(:,1)));
strokepx_in_hMSG=sum(ismember(stroke_map(:,1),stroke_hMSG(:,1)));
strokepx_in_post=sum(ismember(stroke_map(:,1),stroke_post(:,1)));
strokepx_in_pica=sum(ismember(stroke_map(:,1),stroke_pica(:,1)));

strokepx_in_regions=strokepx_in_ant+strokepx_in_vMSG+strokepx_in_hMSG+strokepx_in_post+strokepx_in_pica;


%_________________________
strokepx_per_region_th=0.1;
%_________________________
strokepx_per_region_av=strokepx_per_region_th*strokepx_in_regions; %Absolute value, wieviele Strokepixxel müssen in einer Region sein, damit diese als betroffen gilt


%Überprüfen ob ein sehr großer Teil des Infarktes in nur einer Region liegt
if strokepx_in_ant>=strokepx_per_region_av
    results(casenumber,1,1)=1;
end
if strokepx_in_vMSG>=strokepx_per_region_av
    results(casenumber,2,1)=1;
end
if strokepx_in_hMSG>=strokepx_per_region_av
    results(casenumber,3,1)=1;
end
if strokepx_in_post>=strokepx_per_region_av
    results(casenumber,4,1)=1;
end
if strokepx_in_pica>=strokepx_per_region_av
    results(casenumber,5,1)=1;
end

%________________________
region_covered_by_stroke_aca=0.07;
region_covered_by_stroke_vMSG=0.01;
region_covered_by_stroke_hMSG=0.01;
region_covered_by_stroke_post=0.09;
region_covered_by_stroke_pica=0.01;
%_______________________
%Überprüfen, ob große Teile eines Gebietes mit Strokepixeln bedeckt sind
if strokepx_in_ant>=region_covered_by_stroke_aca*size_ant
    results(casenumber,1,1)=1;
end
if strokepx_in_vMSG>=region_covered_by_stroke_vMSG*size_vMSG
    results(casenumber,2,1)=1;
end
if strokepx_in_hMSG>=region_covered_by_stroke_hMSG*size_hMSG
    results(casenumber,3,1)=1;
end
if strokepx_in_post>=region_covered_by_stroke_post*size_post
    results(casenumber,4,1)=1;
end
if strokepx_in_pica>=region_covered_by_stroke_pica*size_pica
    results(casenumber,5,1)=1;
end
end
save('V:\AG_Neurophotonik\Projekte\Projektion\MATLAB_Programme\results.mat', 'results')


%% Auswertung der berechneten Target Areas

RN=0; % Richtig negativ count
FP=0; % Falsch positiv count
FN=0; % Falsch negativ count
RP=0; % Richtig positiv count
cases_counted=0; % Anzahl der ausgewerteten Cases

%_______________
delta_ratio_threshold=0.5;
%________________
for i=1:numbcases
    

    
    if any(thalamus_cases(:)==i)
        
    else 
    
    if results(i,6,1)>0    
    if results(i,6,1)<delta_ratio_threshold
        cases_counted=cases_counted+1;
         for column=1:5
            if results(i,column,1)==0
                if results(i,column,2)==0
                    RN=RN+1;
                end
                if results(i,column,2)==1
                    FN=FN+1;
                end
            end                
            if results(i,column,1)==1
                if results(i,column,2)==0
                    FP=FP+1;
                end
                if results(i,column,2)==1
                    RP=RP+1;
                end
            end
            
        end
    end
    end
    end
    
end

Youden=(RP/(RP+FN))+(RN/(FP+RN))-1;
sensitivity=((RP)/(RP+FN)); %sens
specificity=((RN)/(RN+FP)); %spec
PPV=RP/(RP+FP); %positive predictive value
NPV=RN/(RN+FN); % negative predictive value
LR=sensitivity/(1-specificity); %positive likelihood ratio
FPR=((FP)/(RN+FP)); %false positive rate

%%Abspeichern der Ergebnisse
ROC_results(count,1)=delta_ratio_threshold;
ROC_results(count,2)=strokepx_per_region_th;
ROC_results(count,3)=region_covered_by_stroke_aca;
ROC_results(count,4)=region_covered_by_stroke_vMSG;
ROC_results(count,5)=region_covered_by_stroke_hMSG;
ROC_results(count,6)=region_covered_by_stroke_post;
ROC_results(count,7)=region_covered_by_stroke_pica;
ROC_results(count,8)=Youden;
ROC_results(count,9)=RP;
ROC_results(count,10)=FP;
ROC_results(count,11)=RN;
ROC_results(count,12)=FN;
ROC_results(count,13)=FPR;
ROC_results(count,14)=sensitivity;
ROC_results(count,15)=specificity;
ROC_results(count,16)=PPV;
ROC_results(count,17)=NPV;
ROC_results(count,18)=LR;
ROC_results(count,19)=cases_counted;
count=count+1;

                       end
                    end
                 end
            end
       end
    end
end


 
% close(findall(0,'Tag','TMWWaitbar'));











