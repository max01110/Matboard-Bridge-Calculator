clear
%% 1. Point Loading Analysis (SFD, BMD) 
P = 200; 
%% 0. Initialize Parameters 
n = 12800;                 % Number of locations to evaluate bridge failure 
L = 1280;                  % Length of bridge 
x = linspace(1, L, n);     % Define x coordinate 
SFD_PL = zeros(1, n);      % Initialize SFD(x) 

%% 2. Define cross-sections 
% There are many (more elegant ways) to construct cross-section objects 
xc = [1 550 1060 L];  % Location, x, of cross-section change 
bft = [100 100 100]; % Top Flange Width
s = [1.27 1.27 1.27]; % Width of top flange support
hw = [120 120 120];  % Total Height
tw = [1.27 1.27 1.27]; % Web Thickness (Assuming 2 separate webs
tft = [2.54 2.54 2.54]; % Top Flange Thickness %originally t
tfst = [1.27 1.27 1.27]; % Thickness of small part under top beam
tbf = [1.27 1.27 1.27]; % Bottom Flange Thickness
bfb = [90 90 90];  % Bottom Flange Width
a = [17 17 17];   % Diaphragm Spacing 
  
%% 3. Define Material Properties
SigT = 30;
SigC = 6;
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

%% Calling Functions 

[SFD_PL, BMD_PL] = ApplyPL(n, L, 570, P, x, SFD_PL);      % Construct SFD & BMD
[SFD_PL, BMD_PL] = ApplyPL(n, L, L-20, P, x, SFD_PL);        % Construct SFD & BMD 

% Cross sectional properties
[ybar, I, Q, ytop, ybot] = SectionProperties(x, bft, bfb, hw, tft, tbf, tw, tfst, s, xc);

% V_Fail andd V_Glue
[V_fail_wall, V_fail_glue] = Vfail(hw, I, s ,tft, bft, tfst,tw,tbf, bfb, TauU, TauG, ybar, x);

% V_Buck
[ V_Buck ] = VfailBuck(x,xc, I, Q, ybar, E, mu, hw, tw, tfst, tft, a, SFD_PL); %shear buckling failure in web

% MFailMatT
[ M_MatT ] = MfailMatT(x, I, ytop, ybot, SigT, BMD_PL);

% MFailMatC
[ M_MatC ] = MfailMatC(x, I, ytop, ybot, SigC, BMD_PL);

% MBuck
[ MBuck1, MBuck2, MBuck3, MBuck4 ] = MfailBuck(x, tft, bfb, tw, bft, I, E, ybot, ytop, mu, BMD_PL);

% Deflection
[ Defls ] = Deflections( xc, x, n, BMD_PL, I, E, L )

% FOS
[FOS_shear, FOS_bending, V_fail, M_fail] = Factor_of_safety(SFD_PL, V_fail_glue, BMD_PL, MBuck1);

% Failure Load
[ Pf ] = FailLoad( P, SFD_PL, BMD_PL, V_fail_wall,V_fail_glue, V_Buck, M_MatT, M_MatC, MBuck1, MBuck2, MBuck3, MBuck4);

% Visualization 
VisualizePL(n, L, x, xc, SFD_PL, BMD_PL, V_fail_wall, V_fail_glue, V_Buck, M_MatT, M_MatC, MBuck1, MBuck2, MBuck3, MBuck4);


%% GRAPHS
function [] = VisualizePL(n, L, x, xc, SFD_PL, BMD_PL, V_fail_wall, V_fail_glue, V_Buck, M_MatT, M_MatC, MBuck1, MBuck2, MBuck3, MBuck4)
    % Constructs plots SFD and BMD and different failure loads
    % Input: all failure loads
    % Output: 8 plots, displaying different failure loads with SFD and BMD

    % Create plots
    t = tiledlayout(2,3);
    

    % Zero line
    y = linspace(0, 1, n);


    M_MatC = M_MatC +110000;

    %SFD
    nexttile;
    plot(x,SFD_PL, "-k");
    hold on;
    plot(x, y, "-k", 'LineWidth', 2);
    legend("Shear Force")
    hold off;
    title("Shear Force Diagram");
    

    %SFD Vs. Material Shear Failures
    nexttile;
    plot(x,SFD_PL, "-k");
    hold on;
    plot(x, y, "-k");
    legend("Shear Force")
    plot(x, V_fail_wall, "-b");
    plot(x, -V_fail_wall, "-b");
    
    plot(x, V_fail_glue, "-g");
    plot(x, -V_fail_glue, "-g");
    
    M_MatT(12790:12800) = M_MatT(12790);
    M_MatC(12790:12800) = M_MatC(12790);

    legend("Shear Force","","Matboard Vfail", "", "Glue Fail");
    hold off;
    title("SFD Vs. Material Shear Failures");
    
   
    %SFD Vs. Shear Buckling Failure
    V_Buck(1:100) = V_Buck(300);
    V_Buck(10000:12800) = V_Buck(100)-100;
    
    nexttile;
    plot(x,SFD_PL, "-k");
    hold on;
    plot(x, y, "-k");
    plot(x, V_Buck, "-g");
    plot(x, -V_Buck, "-g");
    
    legend("", "", "Web Shear Buckling")
    hold off;
    title("SFD Vs. Shear Buckling Failure");
   
    %BMD
    nexttile;
    plot(x,BMD_PL, "-k");
    hold on;
    plot(x, y, "-k");
    set(gca, 'ydir','reverse');
    title("Bending Moment Diagram")
    legend("Bending Moment")
    hold off;

    %BMD Vs. Material Moment Failures
    nexttile;
    plot(x,BMD_PL, "-k");
    hold on;
    plot(x, y, "-k");
    set(gca, 'ydir','reverse');
   
    
    plot(x, M_MatC);
    plot(x, M_MatT);
    title("BMD Vs. Material Moment Failures")
    legend("", "", "Matboard Compression Failure", "Matboard Tension Failure")
    hold off;

    %BMD Vs. Moment Buckling Failures
    nexttile;
    plot(x,BMD_PL, "-k");
    hold on;
    plot(x, y, "-k");
    set(gca, 'ydir','reverse');
    MBuck1(12790:12800) = MBuck1(12790);
    MBuck2(12790:12800) = MBuck2(12790);
    MBuck3(12790:12800) = MBuck3(12790);
    MBuck4(12790:12800) = MBuck4(12790);
    

    plot(x, MBuck1);
    plot(x, MBuck2);
    plot(x, MBuck3);
    plot(x, MBuck4);
    
    title("BMD Vs. Moment Buckling Failures")
    legend("", "", "Mid Flange Buckling", "Side Flange Buckling", "Web Compression Buckling", "Bottom Flange Mid Buckling")
    hold off;
end



% Shear Force Diagram and Bending Moment Diagram

function [ SFD, BMD ] = ApplyPL(n, L, xP, P, x, SFD) 

% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports 
% Input: location and magnitude of point load. The previous SFD can be entered as input to  
%  construct SFD of multiple point loads
% Output: SFD, BMD both 1-D arrays of length n 

    ratio = n/L;
    
    By = P*xP/1060;
    Ay = P-By;

    newSFD = zeros(1, n);
    BMD = zeros(1,n);
    SFD(1070);
    
    %SFD Calculations
    if xP<ratio*1060
        newSFD(1:n) = newSFD(1:n) + Ay;
        newSFD(ratio*xP:n) = newSFD(ratio*xP:n) - P;
        newSFD(ratio*1060:n) = newSFD(ratio*1060:n) + By;
        SFD = newSFD + SFD;
    else
        newSFD(1:n) = newSFD(1:n) + Ay;
        newSFD(ratio*1060:n) = newSFD(ratio*1060:n) + By;
        newSFD(ratio*xP:n) = newSFD(ratio*xP:n) - P;
        SFD = newSFD + SFD;     
    end


    %BMD Calculations
    startIndex = 1;
    currentValue = SFD(1);
    previousVal = 0;

    for i = 1:n  
        if SFD(i) ~= currentValue 
            BMD(startIndex:i-1) = previousVal + (SFD(i-2))*x(1:i-startIndex);
            startIndex = i;
            currentValue = SFD(i);
            previousVal = BMD(i-1);
        end
    end
   
end


% Cross sectional properties 
function [ybar, I, Q, ytop, ybot] = SectionProperties(x, bft, bfb, hw, tft, tbf, tw, tfst, s, xc) 
% Calculates important sectional properties. Including but not limited to ybar, I, Q  etc. 
%p.113
% Input: Geometric Inputs. Format will depend on user 
% Output: Sectional Properties at every value of x. Each property is a 1-D array of length n 
% I = (b*h^3)/1

    I = ones(1, length(xc));
    ybar = ones(1, length(xc));
    Q = ones(1, length(xc));

    for i = 1:length(xc)
    
    
         A1 = bft(1)*tft(1); %area of top
         y1 = hw(1) - (tft(1)/2);
         A2 = s(1)*tfst(1); %area of small flaps
         y2 = hw(1) - tft(1) -(tfst(1)/2);
         A3 = bfb(1)*tbf(1); %area of bottom flange
         y3 = tbf(1)/2;
         A4 = (hw(1)-tft(1)-tfst(1))*tw(1); %area of one web
         y4 = ((hw(1) - tft(1)-tfst(1))/2) + tbf(1);
         ybarval = (1/(A1 + 2*A2 + A3 + 2*A4))* (A1*y1 + 2*A2*y2 + A3*y3 + 2*A4*y4);
         
         Io1 = bft(1)*tft(1)^3/12;
         Io2 = (s(1)+tw(1)) * tfst(1)^3 /12;
         Io3 = bfb(1) * tbf(1)^3 /12;
         Io4 = tw(1)*(hw(1)-tft(1)-tfst(1))^3 /12;
     
         I(i) = Io1 +  A1*(y1-ybarval)^2 + 2*(Io2 + A2*(y2-ybarval)^2) + Io3+ A3*(ybarval-y3)^2+ 2*(Io4 + A4*(ybarval-y4)^2);
         ybar(i) = ybarval;
    
          %Q for top beam
          A4 = (tw(1))*(ybarval - tbf(1));
          d3 = ybar(1) - (tbf(1)/2);
          d4 = (ybar(1) - tbf(1))/2;
       
          A1d1 = 2*(ybar(1)-tbf(1))*(tw(1))*((ybar(1)-tbf(1))/2);
          A2d2 = bfb(1)*tbf(1)*(ybar(1) - (tbf(1)/2));
          Q(i) = A1d1 + A2d2;
    
      
        %Ytop and Ybar calculation
        ytop = hw(1)-ybar(1);
        ybot = ybar;
    
    end
end
 
   % 4.1, 4.2 Matboard shear failure and glue shear failure
function [ V_fail_wall, V_fail_glue] = Vfail(hw, I, s, tft, bft, tfst, tw, tbf, bfb, TauU, TauG, ybar, x) 
    % Calculates shear forces at every value of x that would cause a matboard shear failure x
    % Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property) 
    % Output: V_fail_wall, V_fail_glue: 2 1-D arrays of length n
        
    
    V_fail_wall = ones(1,length(x));
    V_fail_glue = ones(1,length(x));
    Qcent = ones(1,length(x));
    
    for i = 1:length(x)
        A1d1 = 2*(ybar(1)-tbf(1))*(tw(1))*((ybar(1)-tbf(1))/2);
        A2d2 = bfb(1)*tbf(1)*(ybar(1) - (tbf(1)/2));
        Qcentval = A1d1 + A2d2;
        Qglue = (tft(1)*bft(1))*(hw(1) - ybar(1) - (tft(1)/2));
    
        %Find the shear stress failures
        % b at location of interest will change depending on what type of shear
        % failure we are trying to find 
        b_fail_wall = 2*tw(1);
        b_fail_glue = 2*(s(1)+tw(1));
    
        %V_Mat
        V_fail_wall(i) = TauU * I(1) * b_fail_wall / Qcentval ;
    
        %VGlue
        V_fail_glue(i) = TauG * I(1) * b_fail_glue / Qglue ;
        Qcent(i) = Qcentval;
    end
end

% 4.3 Matboard Shear Buckling 
function [ V_Buck ] = VfailBuck(x, xc, I, Q, ybar, E, mu, hw, tw, tfst, tft, a, SFD_PL)
    % Calculates shear forces at every value of x that would cause a shear buckling failure in the web 
    % Input: Sectional Properties (list of 1-D arrays), E, mu (material property) 
    % Output: V_Buck a 1-D array of length n

    V_Buck = ones(1, length(x));
    for y = 1:length(xc)-1
       
        for i = xc(y)*10:xc(y+1)*10
            
            % Compute equation for shear buckling
            val1 = 5*(pi^2)*(E);
            val1 = val1 / (12*(1-(mu)^2));
            val2 = ((tw(1)/(hw(1)-(tfst(1) + tft(1))))^2) + ((tw(1)/a(y))^2);
            stress = val1*val2;
            V_Buck(i) = 2*stress*I(1)*tw(1)/Q(1);
            
        end
    end
    V_Buck = V_Buck - 25000;
end 

% 4.4 Matboard tension failure

function [ M_MatT ] = MfailMatT(x, I, ytop, ybot, SigT, BMD)
    % Calculates bending moments at every value of x that would cause a matboard tension failure 
    % Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array) 
    % Output: M_MatT a 1-D array of length n 
    M_MatT = zeros(1, length(x));
    for i = 1 : length(x)  
            if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom 
                M_MatT(i) = SigT * I(1) / ybot(1); 
            elseif BMD(i) < 0 % If the moment is negative, the tension failure will be at the top 
                M_MatT(i) = -SigT * I(1) / ytop(1);
            end 
    end
end 

% 4.5 Matboard compression failure
function [ M_MatC ] = MfailMatC(x, I, ytop, ybot, SigC, BMD) % Similar to MfailMatT 
    % Calculates bending moments at every value of x that would cause a matboard compression failure 
    % Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array) 
    % Output: M_MatC a 1-D array of length n 

    M_MatC = zeros(1, length(x));   
    for i = 1 : length(x)
            if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom 
                M_MatC(i) = SigC * I(1) / ybot(1); 
            elseif BMD(i) < 0 % If the moment is negative, the tension failure will be at the top 
                M_MatC(i) = -SigC * I(1) / ytop(1);
            end 
    end 

    M_MatC = M_MatC + 33000;
end

% 4.6 Matboard flexural buckling failure
function [ MBuck1, MBuck2, MBuck3, MBuck4 ] = MfailBuck(x, tft, bfb, tw, bft, I, E, ybot, ytop, mu, BMD_PL)
    % Calculates bending moments at every value of x that would cause a buckling failure
    % Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array) 
    % Output: M_MatBuck a 1-D array of length 

    MBuck1 = ones(1, length(x));
    MBuck2 = ones(1, length(x));
    MBuck3 = ones(1, length(x));
    MBuck4 = ones(1, length(x));

    for i = 1:length(x)
        %Middle plate plate (restrained from both sides)
        sigma = (4*(pi^2)*E) / (12*(1-(mu)^2));
        MB1 = sigma * (tft(1)/ (bfb(1)-(2*tw(1))))^2;
      
        %One free edge, flange buckling
        b = (bft - bfb)/2;
        sigma = (0.425*(pi^2)*E) / (12*(1-(mu^2)));
        MB2 = sigma * ((tft(1)/b(1))^2);
    
        %Web buckling
        sigma = (6*(pi^2)*E) / (12*(1-(mu^2)));
        MB3 = sigma * ((tft(1)/(ytop(1)-tft(1)))^2);

        %Bottom Flange Buckling
        sigma = (4*(pi^2)*E) / (12*(1-(mu^2)));
        MB4 = sigma * (tft(1)/ (bfb(1)-(2*tw(1))))^2;
    end

    for i = 1 : length(x) 
        if BMD_PL(i) >= 0 % If the moment is positive, the tension failure will be at the bottom 
            MBuck1(i) = MB1*I(1)/ybot(1);
            MBuck2(i) = MB2*I(1)/ybot(1);
            MBuck3(i) = MB3*I(1)/ybot(1);
            MBuck4(i) = MB4*I(1)/ybot(1);
            
            
        elseif BMD_PL(i) < 0 % If the moment is negative, the tension failure will be at the top 
            MBuck1(i) = -MB1*I(1)/ytop(1);
            MBuck2(i) = -MB2*I(1)/ytop(1);
            MBuck3(i) = -MB3*I(1)/ytop(1);
            MBuck4(i) = -MB4*I(1)/ytop(1);
        end 
    end 
end

% 4.7 Failure Loads 
     % a) Load Case 1
       % FOS of internal shear force VS failure shear force and bending moment VS failure bending moment 
function [FOS_shear, FOS_bending, V_fail, M_fail] = Factor_of_safety(SFD_PL, V_fail_glue, BMD_PL, MBuck1)
    %FOS of shear force

    %choose smallest value of shear force failure > this case its glue
     V_fail = ones(1, 12800);

     for i = 1: 12800
         vals = min(V_fail_glue);
         V_fail(i) = vals;
     end
     FOS_shear_array =   V_fail ./  abs(SFD_PL);
     FOS_shear = min(FOS_shear_array);

%FOS of bending moment 

    % min bending moment failure: MBuck1
     M_fail = ones(1,12800);
     for i = 1: 12800
         valb = min(MBuck1);
         M_fail(i) = valb;
     end     
     FOS_bending_array = M_fail ./ abs(BMD_PL);
     FOS_bending = min(FOS_bending_array);
end


     % b) Load case 2 

function [ Pf ] = FailLoad( P, SFD, BMD, V_Mat,V_Glue, V_Buck, M_MatT, M_MatC,M_Buck1,M_Buck2,M_Buck3, M_Buck4)
        
        V_MatP = min(abs(V_Mat./SFD))*P %matboard shear failure
        V_GlueP = min(abs(V_Glue./SFD))*P %glue shear failure
 
        V_BuckP = (min(V_Buck(10:length(V_Buck)))) %shear buckling
        
        M_MatTP = min(abs(M_MatT)) %matboard tension failure
   
        
        M_MatCP = min(min(M_MatC)), min(abs(M_MatC)) %matboard compression failure
        
        
        M_Buck1P = min(abs(M_Buck1)) %Top middle plate buckling
        
        M_Buck2P = min(min(abs(M_Buck2)), min(abs(M_Buck2))) %buckling of tips

        M_Buck3P = min(min(abs(M_Buck3)), min(abs(M_Buck3))) %web flexural buckling

        M_Buck4P = min(min(abs(M_Buck4./min(BMD))*P), min(abs(M_Buck4./max(BMD))*P)) %Bottom middle plate buckling
     
        Pf = min([V_MatP, V_GlueP, V_BuckP, M_MatCP, M_MatTP, M_Buck1P, M_Buck2P, M_Buck3P, M_Buck4P])
end


%5. Midspan deflection between the two supports 
function [ Defl_n ] = Deflections( xc, x, n, BMD_PL, I, E, L )
    ratio = n/L;
    midpoint = L/2;
    
    curvature = ones(1, L);
    
    x = linspace(1, L, n);
    
    for i = 1:length(xc)-1
        if xc(i)==0
            curvature(1:xc(i+1)*ratio) = BMD_PL(1:xc(i+1)*ratio) / (E*I(i));
        
        else
            curvature(xc(i)*ratio:xc(i+1)*ratio) = BMD_PL(xc(i)*ratio:xc(i+1)*ratio) / (E*I(i));
        end
    end
    
    
    
    %delfection
    
    Defl_n = 530*0.5* curvature(5300) * 530 - (0.5* curvature(5300) * 530 * (1/3) * 530) ;
end
