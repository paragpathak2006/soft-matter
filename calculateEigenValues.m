%% Function that takes in a trial (ld,k1) for given Db and outputs the K matrix
function [L,K] = calculateEigenValues(ld,Db,k1,theta)

h0 = 1;
a = 1;    [~ , ~ ,  ~ , ca ,  ~ ,~ ,~] = matProp(a);    ha = ca*h0/ld;
b = 2;    [~ , ~ ,  ~ , cb ,  ~ ,~ ,~] = matProp(b);    hb = cb*h0/ld;

VaNew = Vfnew(a,ld,Db,k1);
VbNew = Vfnew(b,ld,Db,k1);

[Wa,Za] = eig(VaNew);   ZaVec = Za*ones([6 1]);
[Wb,Zb] = eig(VbNew);   ZbVec = Zb*ones([6 1]);

Qa = Qf(a,ld,Db,k1);
Qb = Qf(b,ld,Db,k1);
% Q = Qa\Qb;
% Qinv = Qb\Qa;
% K = (Qa\Qb)*Wb*expm(Zb*hb)*(Qb\Qa)*Wa*expm(Za*ha)*Wa^-1;
% K1 = (Qa\Qb)*Wb*expm(Zb*hb)*(Qb\Qa)*Wa*expm(Za*ha)*Wa^-1;
% K = (Qa\Qb)*Wb*expm(Zb*hb)*(Qb\Qa)*Wa*expm(Za*ha)*Wa^-1;

    
%          disp('----   -----  ----   ----   -----  -------  -----  ---- ----   -----  ----   ----   -----  -------  -----  ----');
%          fprintf('Input:- lamda = %f , Db = %f, k1 = %f\n',ld,Db,k1);
%          disp('----   -----  ----   ----   -----  -------  -----  ----
%          ----   -----  ----   ----   -----  -------  -----  ----');??

%           Ga = Qa*Wa;
%           Gb = Qb*Wb;
%            K = (Ga\Gb)*expm(Zb*hb)*(Gb\Ga)*expm(Za*ha);
%

Ka = Wb\(Qb\Qa)*Wa*expm(Za*ha);
Kb = Wa\(Qa\Qb)*Wb*expm(Zb*hb);

%     disp('det(VaNew) = ');     disp(det(VaNew))
%     disp('det(VbNew) = ');     disp(det(VbNew))
%
%      disp('expm(Za*ha) = ');     disp(expm(Za*ha));
%      disp('expm(Zb*hb)= ');     disp(expm(Zb*hb));
%
%     disp('Wa^-1 = ');     disp(Wa^-1)
%     disp('Wb^-1 = ');     disp(Wb^-1)
%     disp('Za*ha = ');     disp(Za*ha)
%     disp('Zb*hb = ');     disp(Zb*hb)
%
%      disp('expm(Za*ha) = ');     disp(expm(Za*ha))
%      disp('expm(Zb*hb) = ');     disp(expm(Zb*hb))

%     dKa = det(Ka);
%     dKb = det(Kb);
%     if abs(dKb) > 1000
%         Kb = 6.25*eye(6);
%         dKb = det(Kb);
%     else
K = (Ka)*(Kb);
% K = (Qb\Qa)*Wa*expm(Za*ha)*(Wa\(Qa\Qb))*Wb*expm(Zb*hb)*Wb^-1;
%     end
%     disp(' (Qb\Qa)*Wa  = ');    disp( (Qb\Qa)*Wa)
%     disp('(Qa\Qb)*Wb = ');     disp((Qa\Qb)*Wb)
%
%
%     disp('Ka = ');    disp(dKa);
%     disp('Kb = ');    disp(dKb);
%     disp('K = ');    disp(det(K));
% if abs(det(K)-1)>10000
    %          Zb
    %         K = eye(6);
% end
%     disp('det(K) = ');    disp(det(K));

%   L = eigenOp(K);
e =eig(K);

% VbQ = Q*VbNew*Qinv;
% VaQ = Qinv*VaNew*Q;
% global figType;
% if figType == 1
%     L = det(expm(VbQ*hb) * expm(VaNew*ha)- eye(6,6));
% elseif figType == 2
% 
%     L = det(expm(VbNew*hb) * expm(VaQ*ha)- eye(6,6));
% else
%     L = det(expm(VbNew*hb + VaQ*ha +  (VaQ*ha * VbNew*hb - VbNew*hb * VaQ*ha )^2)- eye(6,6));
%     
% end
L = eigenMatrixDet(K,theta);
%       L1 = abs(e.');%-1.0*ones(6,1)';
%       L = L1(3);
%       L2 = angle(e.');


end
%%
function L = eigenMatrixDet(K,theta)
L = det(K- exp(1i*theta*pi/180)*eye(6,6));
%    if any(imag(L))
%        th = 0;
%    end
end
%%
function Q = Qf(a,ld,Bb,k1)
[B2,B02] = BbTransform(Bb,ld);
[Gm , msmu0 ,  xi , ~ , mu , mu0 ,model] = matProp(a);
[B112,C1221,C1122,A11,A22,B121,B211,B222,C1111,C2121,C1212,C2222,p] = TensorModel(a,msmu0,xi,Gm,ld,B02, mu ,mu0,model);
Q = zeros(6);
Q(1,1) = 1;    
Q(2,3) = 1;      
Q(5,5) = 1;                         
Q(4,6) = -1;Q(4,1) = 1i*k1*(C1122 - C2222 - p); Q(4,5) = B222;
Q(3,2) = C1212;                 Q(3,3) = 1i*k1*(C1221 +  p);        Q(3,4) = B121;                                
Q(6,2) = B121;                  Q(6,3) = 1i*k1*B121;                Q(6,4) = A11;
end

%%
function [V] = Vfnew(a,ld,Bb,k1)
[B2,B02] = BbTransform(Bb,ld);
[Gm , msmu0 ,  xi , ~ , mu , mu0 ,model] = matProp(a);
[B112,C1221,C1122,A11,A22,B121,~,B222,C1111,C2121,C1212,C2222,~] = TensorModel(a,msmu0,xi,Gm,ld, B02, mu ,mu0,model);
B = zeros(6);
B(1,2) = -1;
B(2,1) = k1^2 *(C1122 + C1221 - C1111 );    B(2,5) = 1i*k1* B112;       B(2,6) = -1i*k1;
B(3,1) = 1i*k1;
B(4,1) = (B112 + B121 - B222 )*k1^2;        B(4,5) = -1i*A22*k1;
B(5,4) = 1i*k1;
B(6,2) = (C1221 + C1122 - C2222 )* k1;      B(6,3) = 1i*k1^2*C2121;     B(6,4) = k1*(B121 - B222);
A = eye(6);
A(2,2) = C1212; A(2,4) = B121;
A(4,2) = B121;  A(4,4) = A11;
A(6,6) = 1i;
V = -1*(A\B);
%     a%     A%     B%     [e,ev] = eig(V);%     e%     ev
end
%%
function [B2,B02] = BbTransform(Bb,ld)
global matProVec;
[Gma , ~ ,  ~ , ~ ,mua, mu0 ,~] = matProp(matProVec(13));

B2 = Bb*sqrt(Gma*mu0); %(87) Gma = 1;mu0 = 1; mua = 1;
B02 = B2*ld; %(44)
end
%%

function [Gm , msmu0 ,  xi , c , mu , mu0 ,model] = matProp(a)
[Gma , Gmb , msmu0a , msmu0b ,  xia , xib , ca , cb , mua , mub , mu0 , model,ref]  = getProp();
switch a
    case 1;     Gm = Gma;     msmu0 = msmu0a;       xi = xia;   c = ca;   mu = mua; %Matrix Phase
    case 2;     Gm = Gmb;     msmu0 = msmu0b;       xi = xib;   c = cb;   mu = mub; %Inclusion Phase
end
%     msmu0 =   1.7;    2.16;   2.35;%     xi =      0.75;   0.99;   0.5;%     ci =      0.1;    0.2;    0.3;
end

function [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = TensorModel(a,msmu0,xi,Gm,ld,B0,mu,mu0,model)
switch model
    case 0,     [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = NeoHookeanTensorsWithGamma(msmu0,xi,Gm,ld,B0,mu,mu0);
    case 1,     [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = NeoHookeanTensors(msmu0,xi,Gm,ld,B0,mu,mu0);
    case 2,     [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = LangevinTensors(msmu0,xi,Gm,ld,B0,mu,mu0);
    case 3,     [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = ModifiedNeoHookeanTensors3(msmu0,xi,Gm,ld,B0,mu,mu0);
    case 4,     [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = ModifiedNeoHookeanTensors4(msmu0,xi,Gm,ld,B0,mu,mu0);
    case 5,     [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p]  = getTensors(a);
    otherwise,  [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = DielectricTensors(msmu0,xi,Gm,ld,B0,mu,mu0);
end

end
%% Function that takes in a material properties and gives out Tensor properties given Db,ld
function [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = NeoHookeanTensors(msmu0,~,Gm,ld,B0,mu,mu0)
if msmu0 == 0,    mu=1;end

M112 =  0.0;        A1221 = 0.0;    A1122 = 0.0;
H11 = 1/(mu0*mu);	H22 = 1/(mu0*mu);		%1.0/ep;		%1.0/ep;
M121 = B0/(ld*mu0*mu);M211 = M121;		M222 = 2*B0/(ld*mu0*mu); %D2/ep;										%D2/ep;	%2.0*D2/ep;

A1111 = Gm*ld^2;   A1212 = Gm/ld^2 + B0^2/(ld^2*mu0*mu);  				% mu/ld^2 + D2^2/ep;%mu*ld^2;
A2121 = Gm*ld^2;   A2222 = Gm/ld^2 + B0^2/(ld^2*mu0*mu);    % mu/ld^2 + D2^2/ep;%mu*ld^2;
p = Gm/ld^2 + B0^2*(1-mu)/(ld^2*mu0*mu);

end

%% Function that takes in a material properties and gives out Tensor properties given Db,ld
function [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = ModifiedNeoHookeanTensors3(msmu0,~,Gm,ld,B0,mu,mu0)

if msmu0 == 0,  mu=1; end
% global ldIni;

% if ldIni == 0    
%     ld0 = ld;
% else           
%     ld0 = ldIni;
% end

A1111 = Gm*ld^2;   A1212 = Gm/ld^2 + B0^2/(ld^2*mu0*mu);  % Variable				% mu/ld^2 + D2^2/ep;%mu*ld^2;
B00 = 0;

M112 =  0.0;            A1221 = 0.0;        A1122 = 0.0;
H11 = 1/(mu0*mu);       H22 = 1/(mu0*mu);		%1.0/ep;		%1.0/ep;
M121 = B00/(ld*mu0*mu); M211 = M121;		M222 = 2*B00/(ld*mu0*mu); %D2/ep;										%D2/ep;	%2.0*D2/ep;

A2121 = Gm*ld^2;   A2222 = Gm/ld^2 + B00^2/(ld^2*mu0*mu);    % mu/ld^2 + D2^2/ep;%mu*ld^2;
p = Gm/ld^2+B00^2*(1-mu)/(ld^2*mu0*mu);

end
%%
function [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = ModifiedNeoHookeanTensors4(msmu0,~,Gm,ld,B0,mu,mu0)

if msmu0 == 0,  mu=1; end
% global ldIni;

% if ldIni == 0    
%     ld0 = ld;
% else           
%     ld0 = ldIni;
% end

A1111 = Gm*ld^2;   A2222 = Gm/ld^2 + B0^2/(ld^2*mu0*mu);  % Variable				% mu/ld^2 + D2^2/ep;%mu*ld^2;
B00 = 0;

M112 =  0.0;            A1221 = 0.0;        A1122 = 0.0;
H11 = 1/(mu0*mu);       H22 = 1/(mu0*mu);		%1.0/ep;		%1.0/ep;
M121 = B00/(ld*mu0*mu); M211 = M121;		M222 = 2*B00/(ld*mu0*mu); %D2/ep;										%D2/ep;	%2.0*D2/ep;

A2121 = Gm*ld^2;   A1212 = Gm/ld^2 + B00^2/(ld^2*mu0*mu);    % mu/ld^2 + D2^2/ep;%mu*ld^2;
p = Gm/ld^2+B00^2*(1-mu)/(ld^2*mu0*mu);

end

%%
function [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = LangevinTensors(msmu0,xi,Gm,ld,B0,~,mu0)
%1.0/ep;1.0/ep;D2/ep;2.0*D2/ep;D2/ep;mu*ld^2;mu*ld^2;mu/ld^2 + D2^2/ep;mu/ld^2 + D2^2/ep;

if msmu0 == 0

    H11 = 1/mu0;			    H22 = 1/mu0;		    
    M121 = B0/(ld*mu0);		    M211 = M121;				 M222 = 2*B0/(ld*mu0); 	   
    A1111 = Gm*ld^2;            A2121 = Gm*ld^2;      
    A1212 = Gm/ld^2 + B0^2/(ld^2*mu0);  				     A2222 = Gm/ld^2 + B0^2/(ld^2*mu0);         
    M112 =  0.0;                A1221 = 0.0;                 A1122 = 0.0;    
    p = Gm/ld^2;

    return

end

q = 3*B0*xi/(msmu0*ld); cschp2 = (csch(q))^2;   cothp = coth(q);
p = Gm/ld^2 + msmu0^2/(3*mu0*xi) - B0*msmu0*cothp/(ld*mu0);

M112 =  0.0;    H11 = 1/mu0 + msmu0^2*ld^2/(3*B0^2*xi*mu0) - msmu0*ld*cothp/(B0*mu0);	%1.0/ep;
A1221 = 0.0;    H22 = 1/mu0 - msmu0^2*ld^2/(3*B0^2*xi*mu0) + 3*xi*cschp2/mu0;		%1.0/ep;
A1122 = 0.0;
M121 = B0/(mu0*ld) + msmu0^2*ld/(3*B0*xi*mu0) - msmu0*cothp/mu0;%D2/ep;
M211 =  M121;	M222 = 2*B0/(mu0*ld) - msmu0*cothp/mu0 + 3*B0*xi*cschp2/(ld*mu0); 	%2.0*D2/ep;									%D2/ep;
A1111 = Gm*ld^2;   A1212 = Gm/ld^2 + B0^2/(mu0*ld^2) + msmu0^2/(3*xi*mu0) - B0*msmu0*cothp/(ld*mu0);  				% mu/ld^2 + D2^2/ep;%mu*ld^2;
A2121 = Gm*ld^2;   A2222 = Gm/ld^2 + B0^2/(mu0*ld^2) - msmu0^2/(3*xi*mu0) + 3*B0^2*xi*cschp2/(ld^2*mu0);    % mu/ld^2 + D2^2/ep; %mu*ld^2;

end
%%



%%
function [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p]  = getTensors(a)

global  TensorVec_a TensorVec_b;

if (a == 1),    TensorVec = TensorVec_a;
else,    TensorVec = TensorVec_b;
end

M112 = TensorVec(1);    A1221 = TensorVec(2) ;   A1122 = TensorVec(3)  ;
H11 = TensorVec(4) ;     H22 = TensorVec(5) ;
M121 = TensorVec(6) ;    M211 = TensorVec(7)  ;    M222 = TensorVec(8)  ;
A1111 = TensorVec(9) ;   A2121 = TensorVec(10)  ;  A1212 = TensorVec(11)  ;   A2222 = TensorVec(12)  ;
p = TensorVec(13) ;

end


%% MAEs NeoHookean Tensors With Gamma
function [M112,A1221,A1122,H11,H22,M121,M211,M222,A1111,A2121,A1212,A2222,p] = NeoHookeanTensorsWithGamma(msmu0,~,Gm,ld,B0,mu,mu0)

global r0 r1 r2;
rr0 = r0;    rr1 = r1;     rr2 = r2;

if msmu0 == 0
    mu=1;       r0 = 0;    r1 = 1;     r2 = 0;
end

ld2 = ld^2;
mu0mu = mu0*mu;

M112 =  0.0;        A1221 = 0.0;    A1122 = 0.0;
H11_frac = r1 + r0/ld2 + r2*ld2;
H22_frac = r1 + r0*ld2 + r2/ld2;

M121_frac = r1 + r2*(ld2 + 1/ld2);
M222_frac = r1 + 2*r2/ld2;

A1212_frac = r1 + r2*(ld2 + 2/ld2);
A2222_frac = r1 + r2*6/ld2;
A2121_frac = r2;

H11 = H11_frac/mu0mu;	
H22 = H22_frac/mu0mu;		%1.0/ep;		%1.0/ep;

M121 = B0*M121_frac/(ld*mu0mu);      
M211 = M121;		
M222 = 2*B0*M222_frac/(ld*mu0mu); %D2/ep;										%D2/ep;	%2.0*D2/ep;

A1111 = Gm*ld2;                            
A2121 = Gm*ld2 + B0^2*A2121_frac/(ld2*mu0mu);   

A1212 = Gm/ld2 + B0^2*A1212_frac/(ld2*mu0mu);  				% mu/ld^2 + D2^2/ep;%mu*ld^2;
A2222 = Gm/ld2 + B0^2*A2222_frac/(ld2*mu0mu);    % mu/ld^2 + D2^2/ep;%mu*ld^2;

p_fac = r1 + 2*r2/ld2 - mu;
p = Gm/ld2 + B0^2*p_fac/(ld2*mu0mu);

r0 = rr0;    r1 = rr1;     r2 = rr2;

end

%% Dielectric
function [B112,C1221,C1122,A11,A22,B121,B211,B222,C1111,C2121,C1212,C2222,p] = DielectricTensors(msmu0,~,Gm,ld,B0,mu,mu0)
global r0 r1 r2;
r00 = r0;    r10 = r1;     r20 = r2;
if msmu0 == 0
    mu=1;       r0 = 0;    r1 = 1;     r2 = 0;
end

ld2 = ld^2;

A11_frac = r1 + r0/ld2 + r2*ld2;
A22_frac = r1 + r0*ld2 + r2/ld2;

B121_frac = r1 + r2*(ld2 + 1/ld2);
B222_frac = r1 + 2*r2/ld2;

C1212_frac = r1 + r2*(ld2 + 2/ld2);
C2222_frac = r1 + r2*6/ld2;

ep = mu*mu0;        D2 = B0/ld;

B112 = 0.0;         C1221 = r2*B0^2/ep;        C1122 = 0.0;
A11 = A11_frac/ep;       A22 = A22_frac/ep;
B121 = D2*B121_frac/ep;       B211 = D2*B121_frac/ep;       B222 = 2.0*D2*B222_frac/ep;      
C1111 = Gm*ld2;         C1212 = Gm/ld2 + D2^2*C1212_frac/ep;
C2121 = Gm*ld2 + r2*B0^2/ep;     C2222 = Gm/ld2 + D2^2*C2222_frac/ep;

%Dnd = B0/sqrt(Gm*ep);   %(77)
%p = Gm*(1+Dnd^2)/ld^2;
p_fac = r1 + 2*r2/ld2;
p = Gm/ld2 + B0^2*p_fac/(ep*ld2);

r0 = r00;    r1 = r10;     r2 = r20;

end
