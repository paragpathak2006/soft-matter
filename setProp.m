% Author : Parag Vijayraj Pathak
function setProp(ca,mub,msmu0b,ref,model,Gmb)
    
    mu0 = 1; Gma = 0.0628318531;  msmu0a = 0;  mua = 1;   xia = 0;   xib = (mub-1)/mub;   cb = 1- ca;     

    Gma = 1;
    
    if ref == 'a',  Mat = 1;
    else,           Mat = 2;  end

    global matProVec;
    matProVec = [Gma , Gmb , msmu0a , msmu0b ,  xia , xib , ca , cb , mua , mub , mu0 , model,Mat];
   
end
