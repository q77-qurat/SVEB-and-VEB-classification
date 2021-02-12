
clc;
clear;
load('outdata.mat');

Indexes_qrs=find(QRS>0.12);

for Detect_qrs=length(Indexes_qrs)
   greater_QRS=QRS(Indexes_qrs);%% qrs values which are greater than 0.12
   Rinterval=RRI(Indexes_qrs);
  
  prev_Rinterval=RRI(Indexes_qrs-1);
  Forw_Rinterval=RRI(Indexes_qrs+1);
  
  
find_pvc_prev=find(0.7*prev_Rinterval > Rinterval);
find_pvc_forw=find(0.2*Forw_Rinterval > Rinterval);


if find_pvc_prev~=0
for j=length(find_pvc_prev);
R_PVC=Rinterval(find_pvc_prev);
R_QRS=QRS(find_pvc_prev);

indexesPVC_pre=Indexes_qrs(find_pvc_prev);

disp 'premature ventriculation Signal pre '

end


elseif find_pvc_forw~=0

for j=length(find_pvc_forw);
R_PVC=Rinterval(find_pvc_forw);
R_QRS=QRS(find_pvc_forw);
indexesPVC_for=Indexes_qrs(find_pvc_forw);

disp 'premature ventriculation Signal for'

end 
end


   
   
 solvingRRI=(RRI(1:length(RRI)));
 solvingRRI=solvingRRI';
solvingqrs=(QRS (1:length(QRS)));
solvingqrs=solvingqrs';
solvingHr=(HRC(1:length(HRC)));
solvingHr=solvingHr';

   
end