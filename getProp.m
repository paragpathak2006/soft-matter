% Author: Parag Vijayraj Pathak
function [Gma , Gmb , msmu0a , msmu0b ,  xia , xib , ca , cb , mua , mub , mu0 , model,ref] = getProp()
    
    global matProVec;
    Gma = matProVec(1); 
    Gmb = matProVec(2); 
    msmu0a = matProVec(3) ; 
    msmu0b = matProVec(4);  
    xia  = matProVec(5); 
    xib  = matProVec(6); 
    ca  = matProVec(7); 
    cb  = matProVec(8); 
    mua  = matProVec(9); 
    mub  = matProVec(10); 
    mu0  = matProVec(11); 
    model = matProVec(12);
    ref = matProVec(13);
    
end

% Shear modulus contrast: 20 -- 100
% Initial Susceptibility of the Magnetoactive layer: 0.6 -- 0.999
% Saturation magnetization:: 0.85 -- 2

% mat = [0.6 10.0 10.0];
% mat16 = [.01 16 2.5];
