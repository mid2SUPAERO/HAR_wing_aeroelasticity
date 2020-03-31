% Load factor and steady pitch rate for equilibrium manoeuvre 
n = 1.0;	q_pr = 0;			

% Trim response of a rigid aircraft – requires aircraft data, flight case and derivative codes
% Setting–up and Solving Equations of Motion for the Rigid Aircraft
 
ARigid = - [Z_eta  Z_alpha;  M_eta  M_alpha]; 
CRigid = [1 ; 0];  DRigid = [Z_q ; M_q];  ERigid = [Z_0 ; M_0];
BRigid = inv(ARigid) * (CRigid * (n * W) + DRigid * q_pr + ERigid); Bdeg = BRigid * 180 / pi;
Trim_Elevator_Rigid = Bdeg(1)
Trim_Incidence_Rigid = Bdeg(2)

% Trim response of an elastic aircraft – requires aircraft data, flight case, flexible mode and derivative codes
% Setting–up and Solving Equations of Motion for the Elastic Aircraft 

AElastic = - [Z_eta  Z_alpha  Z_e;  M_eta  M_alpha  M_e;  Q_eta  Q_alpha  Q_e - k_e];  
CElastic = [1 ; 0 ; 0];  DElastic = [Z_q ; M_q ; Q_q];  EElastic = [Z_0 ; M_0 ; Q_0];    
BElastic = inv(AElastic) * (CElastic * (n * W) + DElastic * q_pr + EElastic); Bdeg = BElastic * 180 / pi;
Trim_Elevator_Elastic = Bdeg(1) 
Trim_Incidence_Elastic = Bdeg(2)

% Note that BElastic(3) yields the generalised coordinate for the fuselage bending mode deformation 
% in the trimmed state and so the absolute deformation may be obtained by multiplying by the normalised mode shape.
