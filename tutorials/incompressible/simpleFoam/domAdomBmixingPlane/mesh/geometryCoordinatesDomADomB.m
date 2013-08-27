%
% This Matlab script was written by:
%     Pirooz Moradnia
%     Department of Applied Mechanics
%     Chalmers University of Technology
%     Gothenburg-Sweden
%
% With the collaboration of:
%     Maryse Page, IREQ, Hydro Quebec, Canada
%     Martin Beaudoin, IREQ, Hydro Quebec, Canada
%
%
close all
clear all
clc
% User defined quantities--------------------------------------------------
totalNumberOfStages             =60;
interfaceType                   =char('ggi');        	% Keywords {'stitch','ggi','overlapGgi'}
interfaceTangentialDistribution	=char('NONuniform');    % Keyword {'uniform'}
cyclicBoundaryType              =char('cyclicGgi');     % Keywords {'cyclic','cyclicGgi'}
interfaceRadialuniformity       =char('ifPossible');    % uniform cell size by the interface, Keyword {'ifPossible'}
domBOffsetAngle                 =0;                     % Offset angle in degrees
if ((strcmp(interfaceType,'stitch')==1) && (domBOffsetAngle~=0))
    interfaceType=char('ggi');
    disp('Attention! Selected interface type "stitch" requires the offset angle to be zero');
    disp('The interface type is automatically set to ggi');
end
% Main radii---------------------------------------------------------------
outletDomB                  =4.375;  
bladeDomBCenterRadius       =4.675;  
interfaceRadius             =4.8375;
bladeDomACenterRadius       =5.0;    
inletDomA                   =5.5;    
% z coordinates
zBottom                     =-0.005;
zTop                        =0.005;
% Blades-------------------------------------------------------------------
%
bladeDomALength           =0.15; 
bladeDomAAngle            =-40; 
bladeDomAThickness        =0.025; 
bladeDomAFrameAllowance   =0.025; 
%
bladeDomBLength           =0.15; 
bladeDomBAngle            =-40;   %40; 
bladeDomBThickness        =0.025; 
bladeDomBFrameAllowance   =0.025; 
% -------------------------------------------------------------------------
nR=[36 22 22 8 6  6  6 8 22 22 48  8 8]; % Number of cells in the radial direction of the domain
FR=[1 1 1 1 1 1 1 1 1 1 1];              % Cell ratio in the radial direction of the domain
%
nTDomBZone=[40 10 10 40];          % Number of cells in the tangential direction of the domB domain
FTDomBZone =[1.0 1.0 1.0 1.0];     % Cell ratio in the tangential direction of the domB domain
%
nTDomAZone=[40 10 10 40];         % Number of cells in the tangential direction of the domA domain
FTDomAZone=[1.0 1.0 1.0 1.0];     % Cell ratio in the tangential direction of the domA domain
%
nBT=[16 28 28 16 16 28 28 16];      % Number of cells in the peripheral direction of the domB frame region
nBR=[8 12];                         % Number of cells in the radial direction of the domB frame region
FBT=[1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];   % Cell ratio in the peripheral direction of the domB frame region
FBR=[2.0 1.0];                           % Cell ratio in the radial direction of the domB frame region
%
nAT=[16 28 28 16 16 28 28 16];       % Number of cells in the peripheral direction of the domA frame region
nAR=[8 12];                          % Number of cells in the radial direction of the domA frame region
FAT=[1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];  % Cell ratio in the peripheral direction of the domA frame region
FAR=[2.0 2.0];                          % Cell ratio in the radial direction of the domA frame region
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Angles-------------------------------------------------------------------
bladePassageAngle       =2*pi/totalNumberOfStages;
thetaWedge              =bladePassageAngle/2;
% Domains------------------------------------------------------------------
bladeDomAFrameRadius    =bladeDomALength/2+bladeDomAFrameAllowance;
bladeDomBFrameRadius    =bladeDomBLength/2+bladeDomBFrameAllowance;
%   Rotor------------------------------------------------------------------
R(1)                    =outletDomB;
R(2)                    =bladeDomBCenterRadius-bladeDomBFrameRadius*sin(pi/3);
R(3)                    =bladeDomBCenterRadius;
R(4)                    =bladeDomBCenterRadius+bladeDomBFrameRadius*sin(pi/3);
X(1,:)                  =R(1)*[sin(-thetaWedge) -cos(pi/3)*bladeDomBFrameRadius/bladeDomBCenterRadius 0 cos(pi/3)*bladeDomBFrameRadius/bladeDomBCenterRadius sin(thetaWedge)];
Y(1,:)                  =sqrt(R(1)^2-X(1,:).*X(1,:));
x(1,:)                  =X(1,:)*cos(domBOffsetAngle*pi/180)-Y(1,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
y(1,:)                  =Y(1,:)*cos(domBOffsetAngle*pi/180)+X(1,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
X(2,:)                  =R(2)*[sin(-thetaWedge) -cos(pi/3)*bladeDomBFrameRadius/R(2) 0 cos(pi/3)*bladeDomBFrameRadius/R(2) sin(thetaWedge)];
Y(2,[1 5])              =sqrt(R(2)^2-X(2,[1 5]).*X(2,[1 5]));	% Pirooz 2012.03.29
Y(2,[2 4])              =[bladeDomBCenterRadius-sin(pi/3)*bladeDomBFrameRadius bladeDomBCenterRadius-sin(pi/3)*bladeDomBFrameRadius];	% Pirooz 2012.03.29
Y(2,3)                  =R(3)-bladeDomBFrameRadius;
x(2,:)                  =X(2,:)*cos(domBOffsetAngle*pi/180)-Y(2,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
y(2,:)                  =Y(2,:)*cos(domBOffsetAngle*pi/180)+X(2,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
X(3,:)                  =R(3)*[sin(-thetaWedge) -bladeDomBFrameRadius/R(3) 0 bladeDomBFrameRadius/R(3) sin(thetaWedge)];
Y(3,[1 5])              =sqrt(R(3)^2-X(3,[1 5]).*X(3,[1 5]));   % Pirooz 2012.03.29
Y(3,[2 4])              =[bladeDomBCenterRadius bladeDomBCenterRadius];	% Pirooz 2012.03.29
Y(3,3)                  =0;
x(3,:)                  =X(3,:)*cos(domBOffsetAngle*pi/180)-Y(3,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
y(3,:)                  =Y(3,:)*cos(domBOffsetAngle*pi/180)+X(3,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
X(4,:)                  =R(4)*[sin(-thetaWedge) -cos(pi/3)*bladeDomBFrameRadius/R(4) 0 cos(pi/3)*bladeDomBFrameRadius/R(4) sin(thetaWedge)];
Y(4,[1 5])              =sqrt(R(4)^2-X(4,[1 5]).*X(4,[1 5]));	% Pirooz 2012.03.29
Y(4,[2 4])              =[bladeDomBCenterRadius+sin(pi/3)*bladeDomBFrameRadius bladeDomBCenterRadius+sin(pi/3)*bladeDomBFrameRadius];	% Pirooz 2012.03.29
Y(4,3)                  =R(3)+bladeDomBFrameRadius;
x(4,:)                  =X(4,:)*cos(domBOffsetAngle*pi/180)-Y(4,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
y(4,:)                  =Y(4,:)*cos(domBOffsetAngle*pi/180)+X(4,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
radiusCounter           =size(R,2)+1;
if ((strcmp(interfaceRadialuniformity,'ifPossible')==1) && ((interfaceRadius-(bladeDomBCenterRadius+bladeDomBFrameRadius))>(1/2)*(interfaceRadius-sqrt((bladeDomBCenterRadius+(bladeDomBFrameRadius*sin(pi/3)))^2+(bladeDomBFrameRadius*cos(pi/3))^2))))
    R(radiusCounter)=(interfaceRadius+bladeDomBCenterRadius+bladeDomBFrameRadius)/2;
    X(radiusCounter,:)=R(radiusCounter)*([sin(-thetaWedge) -cos(pi/3)*bladeDomBFrameRadius/bladeDomBCenterRadius 0 cos(pi/3)*bladeDomBFrameRadius/bladeDomBCenterRadius sin(thetaWedge)]);
    Y(radiusCounter,:)=sqrt(R(radiusCounter)^2-X(radiusCounter,:).*X(radiusCounter,:));
    x(radiusCounter,:)=X(radiusCounter,:)*cos(domBOffsetAngle*pi/180)-Y(radiusCounter,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    y(radiusCounter,:)=Y(radiusCounter,:)*cos(domBOffsetAngle*pi/180)+X(radiusCounter,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    radiusCounter=radiusCounter+1;
end
R(radiusCounter)=interfaceRadius;
if ~(strcmp(interfaceTangentialDistribution,'uniform')==1)
    X(radiusCounter,:)=R(radiusCounter)*([sin(-thetaWedge) -cos(pi/3)*bladeDomBFrameRadius/bladeDomBCenterRadius 0 cos(pi/3)*bladeDomBFrameRadius/bladeDomBCenterRadius sin(thetaWedge)]);
    Y(radiusCounter,:)=sqrt(R(radiusCounter)^2-X(radiusCounter,:).*X(radiusCounter,:));
    x(radiusCounter,:)=X(radiusCounter,:)*cos(domBOffsetAngle*pi/180)-Y(radiusCounter,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    y(radiusCounter,:)=Y(radiusCounter,:)*cos(domBOffsetAngle*pi/180)+X(radiusCounter,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
else
    X(radiusCounter,:)=R(radiusCounter)*([sin(-thetaWedge) x(2,2)/R(2) 0 x(2,4)/R(2) sin(thetaWedge)]); % Pirooz 2012.03.29
    Y(radiusCounter,:)=sqrt(R(radiusCounter)^2-X(radiusCounter,:).*X(radiusCounter,:));
    x(radiusCounter,:)=X(radiusCounter,:)*cos(domBOffsetAngle*pi/180)-Y(radiusCounter,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    y(radiusCounter,:)=Y(radiusCounter,:)*cos(domBOffsetAngle*pi/180)+X(radiusCounter,:)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
end
    % -
    % aligned (5,2) and (5,4) nodes on a radial line which the corresponding nodes on the interface
    % -
    counterPT=radiusCounter-1;
    radiusPT=R(counterPT);
    X(counterPT,2)=radiusPT * X(counterPT+1,2)/interfaceRadius;
    Y(counterPT,2)=sqrt( radiusPT*radiusPT - X(counterPT,2)*X(counterPT,2) );
    x(counterPT,2)=X(counterPT,2)*cos(domBOffsetAngle*pi/180)-Y(counterPT,2)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    y(counterPT,2)=Y(counterPT,2)*cos(domBOffsetAngle*pi/180)+X(counterPT,2)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    %
    X(counterPT,4)=radiusPT * X(counterPT+1,4)/interfaceRadius;
    Y(counterPT,4)=sqrt( radiusPT*radiusPT - X(counterPT,4)*X(counterPT,4) );
    x(counterPT,4)=X(counterPT,4)*cos(domBOffsetAngle*pi/180)-Y(counterPT,4)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    y(counterPT,4)=Y(counterPT,4)*cos(domBOffsetAngle*pi/180)+X(counterPT,4)*sin(domBOffsetAngle*pi/180);   % domB domain rotaion
    %
%   Stator-----------------------------------------------------------------
if ~((strcmp(interfaceType,'stitch')==1))
    radiusCounter=radiusCounter+1;
    R(radiusCounter)=interfaceRadius;
    if ~(strcmp(interfaceTangentialDistribution,'uniform')==1)
        x(radiusCounter,:)=R(radiusCounter)*([sin(-thetaWedge) -cos(pi/3)*(bladeDomBFrameRadius+bladeDomAFrameRadius)*0.5/bladeDomACenterRadius 0 cos(pi/3)*(bladeDomBFrameRadius+bladeDomAFrameRadius)*0.5/bladeDomACenterRadius sin(thetaWedge)]);% Pirooz 2012.03.29
    else
        x(radiusCounter,:)=R(radiusCounter)*([sin(-thetaWedge) x(2,2)/R(2) 0 x(2,4)/R(2) sin(thetaWedge)]);% Pirooz 2012.03.29
    end
    y(radiusCounter,:)=sqrt(R(radiusCounter)^2-x(radiusCounter,:).*x(radiusCounter,:));
end
radiusCounter=radiusCounter+1;
if ((strcmp(interfaceRadialuniformity,'ifPossible')==1) && (((bladeDomACenterRadius-bladeDomAFrameRadius)-interfaceRadius)>(1/2)*(sqrt((bladeDomACenterRadius-(bladeDomAFrameRadius*sin(pi/3)))^2+(bladeDomAFrameRadius*cos(pi/3))^2)-interfaceRadius)))
    R(radiusCounter)=(interfaceRadius+bladeDomACenterRadius-bladeDomAFrameRadius)/2;
    x(radiusCounter,:)=R(radiusCounter)*([sin(-thetaWedge) -cos(pi/3)*bladeDomAFrameRadius/bladeDomACenterRadius 0 cos(pi/3)*bladeDomAFrameRadius/bladeDomACenterRadius sin(thetaWedge)]);
    y(radiusCounter,:)=sqrt(R(radiusCounter)^2-x(radiusCounter,:).*x(radiusCounter,:));
    radiusCounter=radiusCounter+1;
end
    % -
    % aligned (8,2) and (8,4) nodes on a radial line which the corresponding nodes on the interface
    % -
    counterPT=radiusCounter-1;
    radiusPT=R(counterPT);
    x(counterPT,2)=radiusPT * x(counterPT-1,2)/interfaceRadius;
    y(counterPT,2)=sqrt( radiusPT*radiusPT - x(counterPT,2)*x(counterPT,2) );
    %
    x(counterPT,4)=radiusPT * x(counterPT-1,4)/interfaceRadius;
    y(counterPT,4)=sqrt( radiusPT*radiusPT - x(counterPT,4)*x(counterPT,4) );
    %
R(radiusCounter)=bladeDomACenterRadius-bladeDomAFrameRadius*sin(pi/3);
x(radiusCounter,:)=R(radiusCounter)*[sin(-thetaWedge) -cos(pi/3)*bladeDomAFrameRadius/R(radiusCounter) 0 cos(pi/3)*bladeDomAFrameRadius/R(radiusCounter) sin(thetaWedge)];
y(radiusCounter,[1 5])=sqrt(R(radiusCounter)^2-x(radiusCounter,[1 5]).*x(radiusCounter,[1 5]));% Pirooz 2012.03.29
y(radiusCounter,[2 4])=[bladeDomACenterRadius-sin(pi/3)*bladeDomAFrameRadius bladeDomACenterRadius-sin(pi/3)*bladeDomAFrameRadius];% Pirooz 2012.03.29
y(radiusCounter,3)=bladeDomACenterRadius-bladeDomAFrameRadius;
radiusCounter=radiusCounter+1;
R(radiusCounter)=bladeDomACenterRadius;
x(radiusCounter,:)=R(radiusCounter)*[sin(-thetaWedge) -bladeDomAFrameRadius/R(radiusCounter) 0 bladeDomAFrameRadius/R(radiusCounter) sin(thetaWedge)];
y(radiusCounter,[1 5])=sqrt(R(radiusCounter)^2-x(radiusCounter,[1 5]).*x(radiusCounter,[1 5]));% Pirooz 2012.03.29
y(radiusCounter,[2 4])=[bladeDomACenterRadius bladeDomACenterRadius];% Pirooz 2012.03.29
y(radiusCounter,3)=0;
radiusCounter=radiusCounter+1;
R(radiusCounter)=bladeDomACenterRadius+bladeDomAFrameRadius*sin(pi/3);
x(radiusCounter,:)=R(radiusCounter)*[sin(-thetaWedge) -cos(pi/3)*bladeDomAFrameRadius/R(radiusCounter) 0 cos(pi/3)*bladeDomAFrameRadius/R(radiusCounter) sin(thetaWedge)];
y(radiusCounter,[1 5])=sqrt(R(radiusCounter)^2-x(radiusCounter,[1 5]).*x(radiusCounter,[1 5]));% Pirooz 2012.03.29
y(radiusCounter,[2 4])=[bladeDomACenterRadius+sin(pi/3)*bladeDomAFrameRadius bladeDomACenterRadius+sin(pi/3)*bladeDomAFrameRadius];% Pirooz 2012.03.29
y(radiusCounter,3)=bladeDomACenterRadius+bladeDomAFrameRadius;
radiusCounter=radiusCounter+1;
R(radiusCounter)=inletDomA;
x(radiusCounter,:)=R(radiusCounter)*([sin(-thetaWedge) -cos(pi/3)*bladeDomAFrameRadius/bladeDomACenterRadius 0 cos(pi/3)*bladeDomAFrameRadius/bladeDomACenterRadius sin(thetaWedge)]);
y(radiusCounter,:)=sqrt(R(radiusCounter)^2-x(radiusCounter,:).*x(radiusCounter,:));
% Number of cells based on the interface type------------------------------
if (strcmp(interfaceTangentialDistribution,'uniform')==1)
        clear nT;
        cellNumberRatio=ceil(10*(cos(pi/3)*(bladeDomBFrameRadius+bladeDomAFrameRadius)*0.5/bladeDomACenterRadius)/(sin(thetaWedge)))/10;
        nT=ceil(sum(nTDomBZone+nTDomAZone)/4*[1-cellNumberRatio cellNumberRatio cellNumberRatio 1-cellNumberRatio]);
        clear nTDomBZone nTDomAZone
        nTDomBZone=nT;
        nTDomAZone=nTDomBZone;
        FTDomBZone=ones(1,4);
        FTDomAZone=FTDomBZone;
else
    if (strcmp(interfaceType,'stitch')==1)
        nTDomBZone=ceil((nTDomBZone+nTDomAZone)/2);
        nTDomAZone=nTDomBZone;
        FTDomBZone=floor((FTDomBZone+FTDomAZone)/2);
        FTDomAZone=FTDomBZone;
    end
end
% Blades-------------------------------------------------------------------
%   Original---------------------------------------------------------------
%       Domain B --------------------------------------------------------------
XB(1,1)=-bladeDomBLength/2;
YB(1,1)=0;
XB(2,1)=-(bladeDomBLength-bladeDomBThickness)/2;
YB(2,1)=-bladeDomBThickness/2;
XB(3,1)=0;
YB(3,1)=-bladeDomBThickness/2;
XB(4,1)=-XB(2,1);
YB(4,1)=YB(2,1);
XB(5,1)=-XB(1,1);
YB(5,1)=0;
XB(6,1)=XB(4,1);
YB(6,1)=-YB(4,1);
XB(7,1)=0;
YB(7,1)=-YB(3,1);
XB(8,1)=XB(2,1);
YB(8,1)=-YB(2,1);
XB(1,3)=-(bladeDomBFrameRadius);
YB(1,3)=0;
XB(2,3)=-(bladeDomBLength)/2;
YB(2,3)=-sqrt(bladeDomBFrameRadius^2-(XB(2,3)^2));
XB(3,3)=0;
YB(3,3)=-bladeDomBFrameRadius;
XB(4,3)=-XB(2,3);
YB(4,3)=YB(2,3);
XB(5,3)=-XB(1,3);
YB(5,3)=0;
XB(6,3)=XB(4,3);
YB(6,3)=-YB(4,3);
XB(7,3)=0;
YB(7,3)=-YB(3,3);
XB(8,3)=XB(2,3);
YB(8,3)=-YB(2,3);
XB(1,2)=-(bladeDomBLength/2+bladeDomBFrameRadius)/2;
YB(1,2)=0;
%XB(2,2)=(XB(2,1)+XB(2,3))/2;
%YB(2,2)=(YB(2,1)+YB(2,3))/2;
XB(2,2)=XB(2,1)+(XB(2,3)-XB(2,1))/3;
YB(2,2)=YB(2,1)+(YB(2,3)-YB(2,1))/3;
XB(3,2)=0;
YB(3,2)=YB(2,2); %-bladeDomBFrameRadius/2;
XB(4,2)=-XB(2,2);
YB(4,2)=YB(2,2);
XB(5,2)=-XB(1,2);
YB(5,2)=0;
XB(6,2)=XB(4,2);
YB(6,2)=-YB(4,2);
XB(7,2)=0;
YB(7,2)=-YB(3,2);
XB(8,2)=XB(2,2);
YB(8,2)=-YB(2,2);
%       Domain A -------------------------------------------------------------
XA(1,1)=-bladeDomALength/2;
YA(1,1)=0;
XA(2,1)=-(bladeDomALength-bladeDomAThickness)/2;
YA(2,1)=-bladeDomAThickness/2;
XA(3,1)=0;
YA(3,1)=-bladeDomAThickness/2;
XA(4,1)=-XA(2,1);
YA(4,1)=YA(2,1);
XA(5,1)=-XA(1,1);
YA(5,1)=0;
XA(6,1)=XA(4,1);
YA(6,1)=-YA(4,1);
XA(7,1)=0;
YA(7,1)=-YA(3,1);
XA(8,1)=XA(2,1);
YA(8,1)=-YA(2,1);
XA(1,3)=-(bladeDomAFrameRadius);
YA(1,3)=0;
XA(2,3)=-(bladeDomALength)/2;
YA(2,3)=-sqrt(bladeDomAFrameRadius^2-(XA(2,3)^2));
XA(3,3)=0;
YA(3,3)=-bladeDomAFrameRadius;
XA(4,3)=-XA(2,3);
YA(4,3)=YA(2,3);
XA(5,3)=-XA(1,3);
YA(5,3)=0;
XA(6,3)=XA(4,3);
YA(6,3)=-YA(4,3);
XA(7,3)=0;
YA(7,3)=-YA(3,3);
XA(8,3)=XA(2,3);
YA(8,3)=-YA(2,3);
XA(1,2)=-(bladeDomALength/2+bladeDomAFrameRadius)/2;
YA(1,2)=0;
%XA(2,2)=(XA(2,1)+XA(2,3))/2;
%YA(2,2)=(YA(2,1)+YA(2,3))/2;
XA(2,2)=XA(2,1)+(XA(2,3)-XA(2,1))/3;
YA(2,2)=YA(2,1)+(YA(2,3)-YA(2,1))/3;
XA(3,2)=0;
YA(3,2)=YA(2,2); %-bladeDomAFrameRadius/2;
XA(4,2)=-XA(2,2);
YA(4,2)=YA(2,2);
XA(5,2)=-XA(1,2);
YA(5,2)=0;
XA(6,2)=XA(4,2);
YA(6,2)=-YA(4,2);
XA(7,2)=0;
YA(7,2)=-YA(3,2);
XA(8,2)=XA(2,2);
YA(8,2)=-YA(2,2);
%   Rotation and translation-----------------------------------------------
%       domB rotation and translation--------------------------------------
XBR(:,:)=XB(:,:)*cos(bladeDomBAngle*pi/180)-YB(:,:)*sin(bladeDomBAngle*pi/180);
YBR(:,:)=YB(:,:)*cos(bladeDomBAngle*pi/180)+XB(:,:)*sin(bladeDomBAngle*pi/180)+bladeDomBCenterRadius;
xR(:,:)=XBR(:,:)*cos(domBOffsetAngle*pi/180)-YBR(:,:)*sin(domBOffsetAngle*pi/180);
yR(:,:)=YBR(:,:)*cos(domBOffsetAngle*pi/180)+XBR(:,:)*sin(domBOffsetAngle*pi/180);
%       domA rotaion and translatio--------------------------------------
xS(:,:)=XA(:,:)*cos(bladeDomAAngle*pi/180)-YA(:,:)*sin(bladeDomAAngle*pi/180);
yS(:,:)=YA(:,:)*cos(bladeDomAAngle*pi/180)+XA(:,:)*sin(bladeDomAAngle*pi/180)+bladeDomACenterRadius;
% Heights------------------------------------------------------------------
height=[zBottom zTop];
heightName=char('B','T');
% calculation loop---------------------------------------------------------

delete blockMeshDict

fid = fopen('blockMeshDict','w');
fprintf(fid,'/*---------------------------------------------------------------------------*\\ \n');
fprintf(fid,'| =========                |                                                 | \n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n');
fprintf(fid,'|  \\    /   O peration     | Version:  1.5                                   | \n');
fprintf(fid,'|   \\  /    A nd           | Web:      http://www.openfoam.org               | \n');
fprintf(fid,'|    \\/     M anipulation  |                                                 | \n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/ \n');
fprintf(fid,' \n');
fprintf(fid,'FoamFile \n');
fprintf(fid,'{ \n');
fprintf(fid,'    version         2.0; \n');
fprintf(fid,'    format          ascii; \n');
fprintf(fid,' \n');
fprintf(fid,'    root            ""; \n');
fprintf(fid,'    case            ""; \n');
fprintf(fid,'    instance        ""; \n');
fprintf(fid,'    local           ""; \n');
fprintf(fid,' \n');
fprintf(fid,'    class           dictionary; \n');
fprintf(fid,'    object          blockMeshDict; \n');
fprintf(fid,'} \n');
fprintf(fid,' \n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //');
fprintf(fid,' \n');
fprintf(fid,'convertToMeters 1; \n');
fprintf(fid,' \n');
fprintf(fid,'vertices \n');
fprintf(fid,'( \n');

vertice=-1;

% vertice------------------------------------------------------------------

for i=1:min((size(nR,2)+1),size(x,1))
    for j=1:size(nTDomBZone,2)+1
        if ~((x(i,j)==0)&&(y(i,j)==0))
            for k=1:size(heightName,1)
                vertice=vertice+1;
                vCount(i,j,k)=vertice;
                fprintf(fid,'\n   (%6.10f	%6.10f	%6.10f)	//	Vertex	R%dT%d%c%	=	%d\n',x(i,j),y(i,j),height(k),i,j,heightName(k),vCount(i,j,k));
            end
        end
    end
end
for iR=1:size(nBT,2)
   for jR=1:size(nBR,2)+1
      for k=1:size(heightName,1)
          vertice=vertice+1;
          vRCount(iR,jR,k)=vertice;
          fprintf(fid,'\n   (%6.10f	%6.10f	%6.10f)	//	Vertex	DomB_R%dT%d%c%	=	%d\n',xR(iR,jR),yR(iR,jR),height(k),iR,jR,heightName(k),vRCount(iR,jR,k));
      end
   end
end
for iS=1:size(nAT,2)
   for jS=1:size(nAR,2)+1
      for k=1:size(heightName,1)
          vertice=vertice+1;
          vSCount(iS,jS,k)=vertice;
          fprintf(fid,'\n   (%6.10f	%6.10f	%6.10f)	//	Vertex	DomA_R%dT%d%c%	=	%d\n',xS(iS,jS),yS(iS,jS),height(k),iS,jS,heightName(k),vSCount(iS,jS,k));
      end
   end
end
fprintf(fid,'\n);\n\n');
fprintf(fid,'////////////////////////////////////////////////////////////////////\n\n');

%Blocks--------------------------------------------------------------------

fprintf(fid,'blocks\n\n');
fprintf(fid,'(\n');
for i=find(R<interfaceRadius)
    for j=1:size(nTDomBZone,2)
        if ~((x(i,j)==0 && y(i,j)==0) || (x(i,j+1)==0 && y(i,j+1)==0) || (x(i+1,j+1)==0 && y(i+1,j+1)==0) || (x(i+1,j)==0 && y(i+1,j)==0))
            fprintf(fid,'	hex (%d %d %d %d %d %d %d %d) domB (%d %d 1) simpleGrading (%6.4f %6.4f 1)\n',vCount(i,j,1),vCount(i,j+1,1),vCount(i+1,j+1,1),vCount(i+1,j,1),vCount(i,j,2),vCount(i,j+1,2),vCount(i+1,j+1,2),vCount(i+1,j,2),nTDomBZone(j),nR(i),FTDomBZone(j),FR(i));
        end
    end
end
for i=[max(find(R==interfaceRadius)) find(R>interfaceRadius & R<inletDomA)]
    for j=1:size(nTDomBZone,2)
        if ~((x(i,j)==0 && y(i,j)==0) || (x(i,j+1)==0 && y(i,j+1)==0) || (x(i+1,j+1)==0 && y(i+1,j+1)==0) || (x(i+1,j)==0 && y(i+1,j)==0))
            fprintf(fid,'	hex (%d %d %d %d %d %d %d %d) domA (%d %d 1) simpleGrading (%6.4f %6.4f 1)\n',vCount(i,j,1),vCount(i,j+1,1),vCount(i+1,j+1,1),vCount(i+1,j,1),vCount(i,j,2),vCount(i,j+1,2),vCount(i+1,j+1,2),vCount(i+1,j,2),nTDomAZone(j),nR(i),FTDomAZone(j),FR(i));
        end
    end
end
for jR=1:size(nBR,2)
    for iR=1:size(nBT,2)-1
        fprintf(fid,'	hex (%d %d %d %d %d %d %d %d) domB (%d %d 1) simpleGrading (%6.4f %6.4f 1)\n',vRCount(iR+1,jR,1),vRCount(iR,jR,1),vRCount(iR,jR+1,1),vRCount(iR+1,jR+1,1),vRCount(iR+1,jR,2),vRCount(iR,jR,2),vRCount(iR,jR+1,2),vRCount(iR+1,jR+1,2),nBT(iR),nBR(jR),FBT(iR),FBR(jR));
    end
    fprintf(fid,'	hex (%d %d %d %d %d %d %d %d) domB (%d %d 1) simpleGrading (%6.4f %6.4f 1)\n',vRCount(1,jR,1),vRCount(8,jR,1),vRCount(8,jR+1,1),vRCount(1,jR+1,1),vRCount(1,jR,2),vRCount(8,jR,2),vRCount(8,jR+1,2),vRCount(1,jR+1,2),nBT(end),nBR(jR),FBT(end),FBR(jR));
end
for jS=1:size(nAR,2)
    for iS=1:size(nAT,2)-1
        fprintf(fid,'	hex (%d %d %d %d %d %d %d %d) domA (%d %d 1) simpleGrading (%6.4f %6.4f 1)\n',vSCount(iS+1,jS,1),vSCount(iS,jS,1),vSCount(iS,jS+1,1),vSCount(iS+1,jS+1,1),vSCount(iS+1,jS,2),vSCount(iS,jS,2),vSCount(iS,jS+1,2),vSCount(iS+1,jS+1,2),nAT(iS),nAR(jS),FAT(iS),FAR(jS));
    end
    fprintf(fid,'	hex (%d %d %d %d %d %d %d %d) domA (%d %d 1) simpleGrading (%6.4f %6.4f 1)\n',vSCount(1,jS,1),vSCount(8,jS,1),vSCount(8,jS+1,1),vSCount(1,jS+1,1),vSCount(1,jS,2),vSCount(8,jS,2),vSCount(8,jS+1,2),vSCount(1,jS+1,2),nAT(end),nAR(jS),FAT(end),FAR(jS));
end
fprintf(fid,');\n\n');
fprintf(fid,'////////////////////////////////////////////////////////////////////\n\n');

%Edges---------------------------------------------------------------------
% Domains------------------------------------------------------------------
fprintf(fid,'edges\n\n');
fprintf(fid,'(\n');
for i=1:min((size(nR,2)+1),size(x,1))
     for j=[1 4]
        xcT(i,j)=0.5*(x(i,j)+x(i,j+1));
        ycT(i,j)=sqrt((R(i))^2-(xcT(i,j))^2);
        for k=1:size(heightName,1)
            if ~(vCount(i,j,k)==0 && vCount(i,j+1,k)==0)
                fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(i,j,k),vCount(i,j+1,k),xcT(i,j),ycT(i,j),height(k));
            end
        end
    end
end
for i=[1 find(R>(bladeDomBCenterRadius+bladeDomBFrameRadius) & R<(bladeDomACenterRadius-bladeDomAFrameRadius)) find(R==inletDomA)]
    for j=[2 3]
        xcT(i,j)=0.5*(x(i,j)+x(i,j+1));
        ycT(i,j)=sqrt((R(i))^2-(xcT(i,j))^2);
        for k=1:size(heightName,1)
            if ~(vCount(i,j,k)==0 && vCount(i,j+1,k)==0)
                fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(i,j,k),vCount(i,j+1,k),xcT(i,j),ycT(i,j),height(k));
            end
        end
    end
end
for j=[2 3]
    XcT(2,j)=0.5*(X(2,j)+X(2,j+1));
    YcT(2,j)=bladeDomBCenterRadius-sqrt((bladeDomBFrameRadius)^2-(XcT(2,j))^2);
    xcT(2,j)=XcT(2,j)*cos(domBOffsetAngle*pi/180)-YcT(2,j)*sin(domBOffsetAngle*pi/180);
    ycT(2,j)=YcT(2,j)*cos(domBOffsetAngle*pi/180)+XcT(2,j)*sin(domBOffsetAngle*pi/180);
    XcT(4,j)=0.5*(X(4,j)+X(4,j+1));
    YcT(4,j)=bladeDomBCenterRadius+sqrt((bladeDomBFrameRadius)^2-(XcT(4,j))^2);
    xcT(4,j)=XcT(4,j)*cos(domBOffsetAngle*pi/180)-YcT(4,j)*sin(domBOffsetAngle*pi/180);
    ycT(4,j)=YcT(4,j)*cos(domBOffsetAngle*pi/180)+XcT(4,j)*sin(domBOffsetAngle*pi/180);
    xcT(size(R,2)-3,j)=0.5*(x(size(R,2)-3,j)+x(size(R,2)-3,j+1));
    ycT(size(R,2)-3,j)=bladeDomACenterRadius-sqrt((bladeDomAFrameRadius)^2-(xcT(size(R,2)-3,j))^2);
    xcT(size(R,2)-1,j)=0.5*(x(size(R,2)-1,j)+x(size(R,2)-1,j+1));
    ycT(size(R,2)-1,j)=bladeDomACenterRadius+sqrt((bladeDomAFrameRadius)^2-(xcT(size(R,2)-1,j))^2);
    for k=1:size(heightName,1)
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(2,j,k),vCount(2,j+1,k),xcT(2,j),ycT(2,j),height(k));
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(4,j,k),vCount(4,j+1,k),xcT(4,j),ycT(4,j),height(k));
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(size(R,2)-3,j,k),vCount(size(R,2)-3,j+1,k),xcT(size(R,2)-3,j),ycT(size(R,2)-3,j),height(k));
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(size(R,2)-1,j,k),vCount(size(R,2)-1,j+1,k),xcT(size(R,2)-1,j),ycT(size(R,2)-1,j),height(k));
    end
end
for j=[2 4]
    XcR(2,j)=0.5*(X(2,j)+X(3,j));
    YcR(2,j)=bladeDomBCenterRadius-sqrt((bladeDomBFrameRadius)^2-(XcR(2,j))^2);
    xcR(2,j)=XcR(2,j)*cos(domBOffsetAngle*pi/180)-YcR(2,j)*sin(domBOffsetAngle*pi/180);
    ycR(2,j)=YcR(2,j)*cos(domBOffsetAngle*pi/180)+XcR(2,j)*sin(domBOffsetAngle*pi/180);
    XcR(3,j)=0.5*(X(3,j)+X(4,j));
    YcR(3,j)=bladeDomBCenterRadius+sqrt((bladeDomBFrameRadius)^2-(XcR(3,j))^2);
    xcR(3,j)=XcR(3,j)*cos(domBOffsetAngle*pi/180)-YcR(3,j)*sin(domBOffsetAngle*pi/180);
    ycR(3,j)=YcR(3,j)*cos(domBOffsetAngle*pi/180)+XcR(3,j)*sin(domBOffsetAngle*pi/180);
    xcR(size(R,2)-2,j)=0.5*(x(size(R,2)-2,j)+x(size(R,2)-1,j));
    ycR(size(R,2)-2,j)=bladeDomACenterRadius+sqrt((bladeDomAFrameRadius)^2-(xcR(size(R,2)-2,j))^2);
    xcR(size(R,2)-3,j)=0.5*(x(size(R,2)-3,j)+x(size(R,2)-2,j));
    ycR(size(R,2)-3,j)=bladeDomACenterRadius-sqrt((bladeDomAFrameRadius)^2-(xcR(size(R,2)-3,j))^2);
    for k=1:size(heightName,1)
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(2,j,k),vCount(3,j,k),xcR(2,j),ycR(2,j),height(k));
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(3,j,k),vCount(4,j,k),xcR(3,j),ycR(3,j),height(k));
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(size(R,2)-2,j,k),vCount(size(R,2)-1,j,k),xcR(size(R,2)-2,j),ycR(size(R,2)-2,j),height(k));
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vCount(size(R,2)-3,j,k),vCount(size(R,2)-2,j,k),xcR(size(R,2)-3,j),ycR(size(R,2)-3,j),height(k));
    end
end
% Blades-------------------------------------------------------------------
%	Domain B ------------------------------------------------------------------
for iR=[1 4 5]
    XBc(iR,1)=0.5*(XB(iR,1)+XB(iR+1,1));
    YBc(iR,1)=sign(YB(iR,1)+YB(iR+1,1))*(bladeDomBThickness/2)*sin(pi/3);
    XBRc(iR,1)=XBc(iR,1)*cos(bladeDomBAngle*pi/180)-YBc(iR,1)*sin(bladeDomBAngle*pi/180);
    YBRc(iR,1)=YBc(iR,1)*cos(bladeDomBAngle*pi/180)+XBc(iR,1)*sin(bladeDomBAngle*pi/180)+bladeDomBCenterRadius;
    xRc(iR,1)=XBRc(iR,1)*cos(domBOffsetAngle*pi/180)-YBRc(iR,1)*sin(domBOffsetAngle*pi/180);
    yRc(iR,1)=YBRc(iR,1)*cos(domBOffsetAngle*pi/180)+XBRc(iR,1)*sin(domBOffsetAngle*pi/180);
    for k=1:size(heightName,1)
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vRCount(iR,1,k),vRCount(iR+1,1,k),xRc(iR,1),yRc(iR,1),height(k));
    end
end
XBc(8,1)=0.5*(XB(8,1)+XB(1,1));
YBc(8,1)=sign(YB(8,1)+YB(1,1))*(bladeDomBThickness/2)*sin(pi/3);
XBRc(8,1)=XBc(8,1)*cos(bladeDomBAngle*pi/180)-YBc(8,1)*sin(bladeDomBAngle*pi/180);
YBRc(8,1)=YBc(8,1)*cos(bladeDomBAngle*pi/180)+XBc(8,1)*sin(bladeDomBAngle*pi/180)+bladeDomBCenterRadius;
xRc(8,1)=XBRc(8,1)*cos(domBOffsetAngle*pi/180)-YBRc(8,1)*sin(domBOffsetAngle*pi/180);
yRc(8,1)=YBRc(8,1)*cos(domBOffsetAngle*pi/180)+XBRc(8,1)*sin(domBOffsetAngle*pi/180);
for k=1:size(heightName,1)
    fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vRCount(8,1,k),vRCount(1,1,k),xRc(8,1),yRc(8,1),height(k));
end
%	Domain A -----------------------------------------------------------------
for iS=[1 4 5]
    XAc(iS,1)=0.5*(XA(iS,1)+XA(iS+1,1));
    YAc(iS,1)=sign(YA(iS,1)+YA(iS+1,1))*(bladeDomAThickness/2)*sin(pi/3);
    xSc(iS,1)=XAc(iS,1)*cos(bladeDomAAngle*pi/180)-YAc(iS,1)*sin(bladeDomAAngle*pi/180);
    ySc(iS,1)=YAc(iS,1)*cos(bladeDomAAngle*pi/180)+XAc(iS,1)*sin(bladeDomAAngle*pi/180)+bladeDomACenterRadius;
    for k=1:size(heightName,1)
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vSCount(iS,1,k),vSCount(iS+1,1,k),xSc(iS,1),ySc(iS,1),height(k));
    end
end
XAc(8,1)=0.5*(XA(8,1)+XA(1,1));
YAc(8,1)=sign(YA(8,1)+YA(1,1))*(bladeDomAThickness/2)*sin(pi/3);
xSc(8,1)=XAc(8,1)*cos(bladeDomAAngle*pi/180)-YAc(8,1)*sin(bladeDomAAngle*pi/180);
ySc(8,1)=YAc(8,1)*cos(bladeDomAAngle*pi/180)+XAc(8,1)*sin(bladeDomAAngle*pi/180)+bladeDomACenterRadius;
for k=1:size(heightName,1)
    fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vSCount(8,1,k),vSCount(1,1,k),xSc(8,1),ySc(8,1),height(k));
end
% Blade frames-------------------------------------------------------------
%	Domain B ------------------------------------------------------------------
 for iR=1:size(nBT,2)-1
     XBc(iR,3)=0.5*(XB(iR,3)+XB(iR+1,3));
     YBc(iR,3)=sign(YB(iR,3)+YB(iR+1,3))*sqrt(bladeDomBFrameRadius^2-(XBc(iR,3))^2);
     XBRc(iR,3)=XBc(iR,3)*cos(bladeDomBAngle*pi/180)-YBc(iR,3)*sin(bladeDomBAngle*pi/180);
     YBRc(iR,3)=YBc(iR,3)*cos(bladeDomBAngle*pi/180)+XBc(iR,3)*sin(bladeDomBAngle*pi/180)+bladeDomBCenterRadius;
     xRc(iR,3)=XBRc(iR,3)*cos(domBOffsetAngle*pi/180)-YBRc(iR,3)*sin(domBOffsetAngle*pi/180);
     yRc(iR,3)=YBRc(iR,3)*cos(domBOffsetAngle*pi/180)+XBRc(iR,3)*sin(domBOffsetAngle*pi/180);
     for k=1:size(heightName,1)
         fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vRCount(iR,3,k),vRCount(iR+1,3,k),xRc(iR,3),yRc(iR,3),height(k));
     end
 end
 XBc(8,3)=0.5*(XB(8,3)+XB(1,3));
 YBc(8,3)=sign(YB(8,3)+YB(1,3))*sqrt(bladeDomBFrameRadius^2-(XBc(8,3))^2);
 XBRc(8,3)=XBc(8,3)*cos(bladeDomBAngle*pi/180)-YBc(8,3)*sin(bladeDomBAngle*pi/180);
 YBRc(8,3)=YBc(8,3)*cos(bladeDomBAngle*pi/180)+XBc(8,3)*sin(bladeDomBAngle*pi/180)+bladeDomBCenterRadius;
 xRc(8,3)=XBRc(8,3)*cos(domBOffsetAngle*pi/180)-YBRc(8,3)*sin(domBOffsetAngle*pi/180);
 yRc(8,3)=YBRc(8,3)*cos(domBOffsetAngle*pi/180)+XBRc(8,3)*sin(domBOffsetAngle*pi/180);
 for k=1:size(heightName,1)
     fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vRCount(8,3,k),vRCount(1,3,k),xRc(8,3),yRc(8,3),height(k));
 end
%	Domain A -----------------------------------------------------------------
for iS=1:size(nAT,2)-1
    XAc(iS,3)=0.5*(XA(iS,3)+XA(iS+1,3));
    YAc(iS,3)=sign(YA(iS,3)+YA(iS+1,3))*sqrt(bladeDomAFrameRadius^2-(XAc(iS,3))^2);
    xSc(iS,3)=XAc(iS,3)*cos(bladeDomAAngle*pi/180)-YAc(iS,3)*sin(bladeDomAAngle*pi/180);
    ySc(iS,3)=YAc(iS,3)*cos(bladeDomAAngle*pi/180)+XAc(iS,3)*sin(bladeDomAAngle*pi/180)+bladeDomACenterRadius;
    for k=1:size(heightName,1)
        fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vSCount(iS,3,k),vSCount(iS+1,3,k),xSc(iS,3),ySc(iS,3),height(k));
    end
end
XAc(8,3)=0.5*(XA(8,3)+XA(1,3));
YAc(8,3)=sign(YA(8,3)+YA(1,3))*sqrt(bladeDomAFrameRadius^2-(XAc(8,3))^2);
xSc(8,3)=XAc(8,3)*cos(bladeDomAAngle*pi/180)-YAc(8,3)*sin(bladeDomAAngle*pi/180);
ySc(8,3)=YAc(8,3)*cos(bladeDomAAngle*pi/180)+XAc(8,3)*sin(bladeDomAAngle*pi/180)+bladeDomACenterRadius;
for k=1:size(heightName,1)
    fprintf(fid,'	arc %d %d (%6.10f %6.10f %6.10f)\n',vSCount(8,3,k),vSCount(1,3,k),xSc(8,3),ySc(8,3),height(k));
end
fprintf(fid,'\n');
fprintf(fid,');\n\n');
fprintf(fid,'////////////////////////////////////////////////////////////////////\n\n');

%Patches-------------------------------------------------------------------

fprintf(fid,'patches\n\n');
fprintf(fid,'(\n');
fprintf(fid,'wall bladeDomB\n');
fprintf(fid,'   (\n');
for iR=1:size(nBT,2)-1
    fprintf(fid,'   (%d %d %d %d)\n',vRCount(iR+1,1,1),vRCount(iR,1,1),vRCount(iR,1,2),vRCount(iR+1,1,2));
end
fprintf(fid,'   (%d %d %d %d)\n',vRCount(1,1,1),vRCount(8,1,1),vRCount(8,1,2),vRCount(1,1,2));
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fprintf(fid,'wall bladeDomA\n');
fprintf(fid,'   (\n');
for iS=1:size(nAT,2)-1
    fprintf(fid,'   (%d %d %d %d)\n',vSCount(iS+1,1,1),vSCount(iS,1,1),vSCount(iS,1,2),vSCount(iS+1,1,2));
end
fprintf(fid,'   (%d %d %d %d)\n',vSCount(1,1,1),vSCount(8,1,1),vSCount(8,1,2),vSCount(1,1,2));
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fprintf(fid,'patch outletDomB\n');
fprintf(fid,'   (\n');
for j=1:size(nTDomBZone,2)
    fprintf(fid,'   (%d %d %d %d)\n',vCount(1,j,1),vCount(1,j+1,1),vCount(1,j+1,2),vCount(1,j,2));
end
fprintf(fid,'   )\n\n');

fprintf(fid,'patch inletDomA\n');
fprintf(fid,'   (\n');
for j=1:size(nTDomBZone,2)
    fprintf(fid,'   (%d %d %d %d)\n',vCount(size(R,2),j+1,1),vCount(size(R,2),j,1),vCount(size(R,2),j,2),vCount(size(R,2),j+1,2));
end
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if ~((strcmp(interfaceType,'stitch')==1) && (domBOffsetAngle==0))    
    fprintf(fid,'%s outletDomA\n',interfaceType);
    fprintf(fid,'   (\n');
    for j=1:size(nTDomBZone,2)
        fprintf(fid,'   (%d %d %d %d)\n',vCount(max(find(R==interfaceRadius)),j,1),vCount(max(find(R==interfaceRadius)),j+1,1),vCount(max(find(R==interfaceRadius)),j+1,2),vCount(max(find(R==interfaceRadius)),j,2));
    end
    fprintf(fid,'   )\n\n');
    %//////////////////////////////////////////////////////////////////////
    %//////////////////////////////////////////////////////////////////////
    fprintf(fid,'%s	inletDomB\n',interfaceType);
    fprintf(fid,'   (\n');
    for j=1:size(nTDomBZone,2)
        fprintf(fid,'   (%d %d %d %d)\n',vCount(min(find(R==interfaceRadius)),j+1,1),vCount(min(find(R==interfaceRadius)),j,1),vCount(min(find(R==interfaceRadius)),j,2),vCount(min(find(R==interfaceRadius)),j+1,2));
    end
    fprintf(fid,'   )\n\n');
end

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fprintf(fid,'wall	ggi1DomB\n');
fprintf(fid,'   (\n');
fprintf(fid,'   (%d %d %d %d)\n',vCount(2,3,1),vCount(2,2,1),vCount(2,2,2),vCount(2,3,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(2,4,1),vCount(2,3,1),vCount(2,3,2),vCount(2,4,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(3,4,1),vCount(2,4,1),vCount(2,4,2),vCount(3,4,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(4,4,1),vCount(3,4,1),vCount(3,4,2),vCount(4,4,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(4,3,1),vCount(4,4,1),vCount(4,4,2),vCount(4,3,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(4,2,1),vCount(4,3,1),vCount(4,3,2),vCount(4,2,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(3,2,1),vCount(4,2,1),vCount(4,2,2),vCount(3,2,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(2,2,1),vCount(3,2,1),vCount(3,2,2),vCount(2,2,2));
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

outletRow=min((size(nR,2)),(size(x,1)-1));
fprintf(fid,'wall	ggi1DomA\n');
fprintf(fid,'   (\n');
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow-2,3,1),vCount(outletRow-2,2,1),vCount(outletRow-2,2,2),vCount(outletRow-2,3,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow-2,4,1),vCount(outletRow-2,3,1),vCount(outletRow-2,3,2),vCount(outletRow-2,4,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow-2,4,1),vCount(outletRow-1,4,1),vCount(outletRow-1,4,2),vCount(outletRow-2,4,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow-1,4,1),vCount(outletRow,4,1),vCount(outletRow,4,2),vCount(outletRow-1,4,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow,3,1),vCount(outletRow,4,1),vCount(outletRow,4,2),vCount(outletRow,3,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow,2,1),vCount(outletRow,3,1),vCount(outletRow,3,2),vCount(outletRow,2,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow-1,2,1),vCount(outletRow-2,2,1),vCount(outletRow-2,2,2),vCount(outletRow-1,2,2));
fprintf(fid,'   (%d %d %d %d)\n',vCount(outletRow,2,1),vCount(outletRow-1,2,1),vCount(outletRow-1,2,2),vCount(outletRow,2,2));
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fprintf(fid,'wall	ggi2DomB\n');
fprintf(fid,'   (\n');
for iR=1:size(nBT,2)-1
    fprintf(fid,'   (%d %d %d %d)\n',vRCount(iR,3,1),vRCount(iR+1,3,1),vRCount(iR+1,3,2),vRCount(iR,3,2));
end
fprintf(fid,'   (%d %d %d %d)\n',vRCount(8,3,1),vRCount(1,3,1),vRCount(1,3,2),vRCount(8,3,2));
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fprintf(fid,'wall	ggi2DomA\n');
fprintf(fid,'   (\n');
for iS=1:size(nAT,2)-1
    fprintf(fid,'   (%d %d %d %d)\n',vSCount(iS,3,1),vSCount(iS+1,3,1),vSCount(iS+1,3,2),vSCount(iS,3,2));
end
fprintf(fid,'   (%d %d %d %d)\n',vSCount(8,3,1),vSCount(1,3,1),vSCount(1,3,2),vSCount(8,3,2));
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fprintf(fid,'empty topAndBottom\n');
fprintf(fid,'   (\n');
for i=[find(R<interfaceRadius) max(find(R==interfaceRadius)) find(R>interfaceRadius & R<inletDomA)]
    for j=1:size(nTDomBZone,2)
        if ~((x(i,j)==0 && y(i,j)==0) || (x(i,j+1)==0 && y(i,j+1)==0) || (x(i+1,j+1)==0 && y(i+1,j+1)==0) || (x(i+1,j)==0 && y(i+1,j)==0))
            fprintf(fid,'   (%d %d %d %d)\n',vCount(i,j,2),vCount(i,j+1,2),vCount(i+1,j+1,2),vCount(i+1,j,2));
            fprintf(fid,'   (%d %d %d %d)\n',vCount(i+1,j,1),vCount(i+1,j+1,1),vCount(i,j+1,1),vCount(i,j,1));
        end
    end    
end
for jR=1:size(nBR,2)
    for iR=1:size(nBT,2)-1
        fprintf(fid,'   (%d %d %d %d)\n',vRCount(iR,jR,1),vRCount(iR+1,jR,1),vRCount(iR+1,jR+1,1),vRCount(iR,jR+1,1));
        fprintf(fid,'   (%d %d %d %d)\n',vRCount(iR+1,jR,2),vRCount(iR,jR,2),vRCount(iR,jR+1,2),vRCount(iR+1,jR+1,2));
    end
    fprintf(fid,'   (%d %d %d %d)\n',vRCount(8,jR,1),vRCount(1,jR,1),vRCount(1,jR+1,1),vRCount(8,jR+1,1));
    fprintf(fid,'   (%d %d %d %d)\n',vRCount(1,jR,2),vRCount(8,jR,2),vRCount(8,jR+1,2),vRCount(1,jR+1,2));
end
for jS=1:size(nAR,2)
    for iS=1:size(nAT,2)-1
        fprintf(fid,'   (%d %d %d %d)\n',vSCount(iS,jS,1),vSCount(iS+1,jS,1),vSCount(iS+1,jS+1,1),vSCount(iS,jS+1,1));
        fprintf(fid,'   (%d %d %d %d)\n',vSCount(iS+1,jS,2),vSCount(iS,jS,2),vSCount(iS,jS+1,2),vSCount(iS+1,jS+1,2));
    end
    fprintf(fid,'   (%d %d %d %d)\n',vSCount(8,jS,1),vSCount(1,jS,1),vSCount(1,jS+1,1),vSCount(8,jS+1,1));
    fprintf(fid,'   (%d %d %d %d)\n',vSCount(1,jS,2),vSCount(8,jS,2),vSCount(8,jS+1,2),vSCount(1,jS+1,2));
end
fprintf(fid,'   )\n\n');

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if strcmp(cyclicBoundaryType,'cyclic')==1
    fprintf(fid,'cyclic perDomB\n');
    fprintf(fid,'   (\n');
    for i=[find(R<interfaceRadius)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i+1,1,1),vCount(i,1,1),vCount(i,1,2),vCount(i+1,1,2));
    end
    for i=[find(R<interfaceRadius)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,2),vCount(i,size(nTDomBZone,2)+1,2));
    end
    fprintf(fid,'   )\n\n');
    fprintf(fid,'cyclic	perDomA\n');
    fprintf(fid,'   (\n');
    for i=[max(find(R==interfaceRadius)) find(R>interfaceRadius & R<inletDomA)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i+1,1,1),vCount(i,1,1),vCount(i,1,2),vCount(i+1,1,2));
    end
    for i=[max(find(R==interfaceRadius)) find(R>interfaceRadius & R<inletDomA)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,2),vCount(i,size(nTDomBZone,2)+1,2));
    end
    fprintf(fid,'   )\n\n');
elseif strcmp(cyclicBoundaryType,'cyclicGgi')==1
    fprintf(fid,'wall	per1DomB\n');
    fprintf(fid,'   (\n');
    for i=[find(R<interfaceRadius)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i+1,1,1),vCount(i,1,1),vCount(i,1,2),vCount(i+1,1,2));        
    end
    fprintf(fid,'   )\n\n');
    fprintf(fid,'wall	per2DomB\n');
    fprintf(fid,'   (\n');
    for i=[find(R<interfaceRadius)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,2),vCount(i,size(nTDomBZone,2)+1,2));
    end
    fprintf(fid,'   )\n\n');
    fprintf(fid,'wall	per1DomA\n');
    fprintf(fid,'   (\n');
    for i=[max(find(R==interfaceRadius)) find(R>interfaceRadius & R<inletDomA)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i+1,1,1),vCount(i,1,1),vCount(i,1,2),vCount(i+1,1,2));
    end
    fprintf(fid,'   )\n\n');
    fprintf(fid,'wall	per2DomA\n');
    fprintf(fid,'   (\n');
    for i=[max(find(R==interfaceRadius)) find(R>interfaceRadius & R<inletDomA)];
        fprintf(fid,'   (%d %d %d %d)\n',vCount(i,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,1),vCount(i+1,size(nTDomBZone,2)+1,2),vCount(i,size(nTDomBZone,2)+1,2));
    end
    fprintf(fid,'   )\n\n');
end

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fprintf(fid,');\n\n');
fprintf(fid,'////////////////////////////////////////////////////////////////////\n\n');
fclose(fid);
