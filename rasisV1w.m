function dSt=rasisV1w(t,St)
%Differential equations 

global B wm wi1 wi2 u um beta b1 phi1 b2 phi2 d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv reintro wA; %immig deaths 


lamda=(1+b1*cos(2*pi*(t-phi1*12)/12)+b2*cos(2*pi*(t-phi2*12)/6))*beta*(St(2*al+1:3*al)+ri2*St(5*al+1:6*al)+ri3*St(8*al+1:9*al)+St(16*al+1:17*al)+ri2*St(19*al+1:20*al)+ri3*St(22*al+1:23*al))/sum(St);

for i=1:al
    dSt(i,1)=log(1+B(round(t),i))/12*sum(St(:))*(1-v1(round(t),1))  - (u(i)+um(i))*St(i) - wm*St(i); %dM/dt  
    
    dSt(i+al,1)=wm*St(i+14*al) - St(i+al)*(sum(lamda(i,:))+reintro) - (u(i)+um(i))*St(i+al); %dS0/dt 
    dSt(i+2*al,1)=St(i+al)*(sum(lamda(i,:))+reintro) - (d1+u(i)+um(i))*St(i+2*al); %dI1/dt
    dSt(i+3*al,1)=d1*St(i+2*al) - (wi1+u(i)+um(i))*St(i+3*al); %dR1/dt
    dSt(i+4*al,1)=wi1*St(i+3*al) - St(i+4*al)*rr1*(sum(lamda(i,:))+reintro) - (u(i)+um(i))*St(i+4*al) + wA*St(i+7*al); %dS1/dt 
    dSt(i+5*al,1)=St(i+4*al)*rr1*(sum(lamda(i,:))+reintro) - (d2+u(i)+um(i))*St(i+5*al); %dI2/dt
    dSt(i+6*al,1)=d2*St(i+5*al) - (wi1+u(i)+um(i))*St(i+6*al); %dR2/dt
    dSt(i+7*al,1)=wi1*St(i+6*al) + wi2*St(i+9*al) - St(i+7*al)*rr2*(sum(lamda(i,:))+reintro) - (u(i)+um(i)+wA)*St(i+7*al); %dSR/dt
    dSt(i+8*al,1)=St(i+7*al)*rr2*(sum(lamda(i,:))+reintro) - (d2+u(i)+um(i))*St(i+8*al); %dAI/dt
    dSt(i+9*al,1)=d2*St(i+8*al) - (wi2+u(i)+um(i))*St(i+9*al); %dR/dt
    
    dSt(i+10*al,1)=wm*St(i) - (u(i)+um(i))*St(i+10*al) - wm*St(i+10*al); %dM2/dt
    dSt(i+11*al,1)=wm*St(i+10*al) - (u(i)+um(i))*St(i+11*al) - wm*St(i+11*al); %dM3/dt
    dSt(i+12*al,1)=wm*St(i+11*al) - (u(i)+um(i))*St(i+12*al) - wm*St(i+12*al); %dM4/dt
    dSt(i+13*al,1)=wm*St(i+12*al) - (u(i)+um(i))*St(i+13*al) - wm*St(i+13*al); %dM5/dt
    dSt(i+14*al,1)=wm*St(i+13*al) - (u(i)+um(i))*St(i+14*al) - wm*St(i+14*al); %dM6/dt
    
    dSt(i+15*al,1)=wm*St(i+29*al) + wv*St(i+30*al) - St(i+15*al)*(sum(lamda(i,:))+reintro) - (u(i)+um(i))*St(i+15*al); %dSV0/dt 
    dSt(i+16*al,1)=St(i+15*al)*(sum(lamda(i,:))+reintro) - (d1+u(i)+um(i))*St(i+16*al); %dIV0/dt
    dSt(i+17*al,1)=d1*St(i+16*al) - (wi1+u(i)+um(i))*St(i+17*al); %dVR0/dt
    dSt(i+18*al,1)=wi1*St(i+17*al) + wv*St(i+31*al) - St(i+18*al)*rr1*(sum(lamda(i,:))+reintro) - (u(i)+um(i))*St(i+18*al); %dSV1/dt 
    dSt(i+19*al,1)=St(i+18*al)*rr1*(sum(lamda(i,:))+reintro) - (d2+u(i)+um(i))*St(i+19*al); %dIV1/dt
    dSt(i+20*al,1)=d2*St(i+19*al) - (wi1+u(i)+um(i))*St(i+20*al); %dVR1/dt
    dSt(i+21*al,1)=wi1*St(i+20*al) + wi2*St(i+23*al) + wv*St(i+32*al) - St(i+21*al)*rr2*(sum(lamda(i,:))+reintro) - (u(i)+um(i))*St(i+21*al); %dSV2/dt 
    dSt(i+22*al,1)=St(i+21*al)*rr2*(sum(lamda(i,:))+reintro) - (d2+u(i)+um(i))*St(i+22*al); %dIV2/dt
    dSt(i+23*al,1)=d2*St(i+22*al) - (wi2+u(i)+um(i))*St(i+23*al); %dVR2/dt
    
    dSt(i+24*al,1)=((1-sc1)*v1(round(t),1)+(1-sc1-sc2)*v2(round(t),1))*log(1+B(round(t),i))/12*sum(St(:)) - (wm+u(i)+um(i))*St(i+24*al); %dMV1/dt
    dSt(i+25*al,1)=wm*St(i+24*al) - (wm+u(i)+um(i))*St(i+25*al); %dMV2/dt
    dSt(i+26*al,1)=wm*St(i+25*al) - (wm+u(i)+um(i))*St(i+26*al); %dMV3/dt
    dSt(i+27*al,1)=wm*St(i+26*al) - (wm+u(i)+um(i))*St(i+27*al); %dMV4/dt
    dSt(i+28*al,1)=wm*St(i+27*al) - (wm+u(i)+um(i))*St(i+28*al); %dMV5/dt
    dSt(i+29*al,1)=wm*St(i+28*al) - (wm+u(i)+um(i))*St(i+29*al); %dMV6/dt
    
    dSt(i+30*al,1)=(sc1*v1(round(t),1)+sc2*(1-sc1)*v2(round(t),1))*log(1+B(round(t),i))/12*sum(St(:)) - (wv+u(i)+um(i))*St(i+30*al); %dV0/dt
    dSt(i+31*al,1)=sc1*sc2*v2(round(t),1)*log(1+B(round(t),i))/12*sum(St(:)) - (wv+u(i)+um(i))*St(i+31*al); %dV1/dt 
    dSt(i+32*al,1)=sc1*sc2*sc3*v3(round(t),1)*log(1+B(round(t),i))/12*sum(St(:)) - (wv+u(i)+um(i))*St(i+32*al); %dV2/dt
    
    if i>1 %aging of individuals from one age group to the next and vaccination
        dSt(i,1)=dSt(i,1) + (1-v1(round(t),i))*u(i-1)*St(i-1);
        dSt(i+al,1)=dSt(i+al,1) + (1-v1(round(t),i))*u(i-1)*St(i+al-1);
        dSt(i+2*al,1)=dSt(i+2*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+2*al-1);
        dSt(i+3*al,1)=dSt(i+3*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+3*al-1);
        dSt(i+4*al,1)=dSt(i+4*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+4*al-1); 
        dSt(i+5*al,1)=dSt(i+5*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+5*al-1);
        dSt(i+6*al,1)=dSt(i+6*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+6*al-1);
        dSt(i+7*al,1)=dSt(i+7*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+7*al-1);
        dSt(i+8*al,1)=dSt(i+8*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+8*al-1);
        dSt(i+9*al,1)=dSt(i+9*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+9*al-1);
                
        dSt(i+10*al,1)=dSt(i+10*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+10*al-1);
        dSt(i+11*al,1)=dSt(i+11*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+11*al-1);
        dSt(i+12*al,1)=dSt(i+12*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+12*al-1);
        dSt(i+13*al,1)=dSt(i+13*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+13*al-1);
        dSt(i+14*al,1)=dSt(i+14*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+14*al-1);
        
        dSt(i+15*al,1)=dSt(i+15*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+15*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*(St(i+al-1)); %SV0
        dSt(i+16*al,1)=dSt(i+16*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+16*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*St(i+2*al-1); %IV0
        dSt(i+17*al,1)=dSt(i+17*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+17*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*St(i+3*al-1); %RV0
        dSt(i+18*al,1)=dSt(i+18*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+18*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*St(i+4*al-1); %SV1
        dSt(i+19*al,1)=dSt(i+19*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+19*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*St(i+5*al-1) + sc1*v1(round(t),i)*u(i-1)*St(i+2*al-1) + (sc2n*v2(round(t),i)+sc3n*v3(round(t),i))*u(i-1)*St(i+16*al-1); %IV1
        dSt(i+20*al,1)=dSt(i+20*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+20*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*St(i+6*al-1); %RV1
        dSt(i+21*al,1)=dSt(i+21*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+21*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*St(i+7*al-1); %SV2
        dSt(i+22*al,1)=dSt(i+22*al,1) + u(i-1)*St(i+22*al-1) + v1(round(t),i)*u(i-1)*St(i+8*al-1) + sc1*v1(round(t),i)*u(i-1)*St(i+5*al-1) + (sc2*v2(round(t),i)+sc3*v3(round(t),i))*u(i-1)*St(i+19*al-1); %IV2
        dSt(i+23*al,1)=dSt(i+23*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+23*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*St(i+9*al-1); %RV2

        dSt(i+24*al,1)=dSt(i+24*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+24*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*(St(i-1)); %MV0
        dSt(i+25*al,1)=dSt(i+25*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+25*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*(St(i+10*al-1)); %
        dSt(i+26*al,1)=dSt(i+26*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+26*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*(St(i+11*al-1)); %
        dSt(i+27*al,1)=dSt(i+27*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+27*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*(St(i+12*al-1)); %
        dSt(i+28*al,1)=dSt(i+28*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+28*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*(St(i+13*al-1)); %
        dSt(i+29*al,1)=dSt(i+29*al,1) + (1-sc2n*v2(round(t),i))*(1-sc3n*v3(round(t),i))*u(i-1)*St(i+29*al-1) + (1-sc1)*v1(round(t),i)*u(i-1)*(St(i+14*al-1)); %
        
        dSt(i+30*al,1)=dSt(i+30*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+30*al-1) + sc1*v1(round(t),i)*u(i-1)*(St(i-1)+St(i+al-1)+St(i+10*al-1)+St(i+11*al-1)+St(i+12*al-1)+St(i+13*al-1)+St(i+14*al-1)) + (sc2n*v2(round(t),i)+sc3n*v3(round(t),i))*u(i-1)*(St(i+15*al-1)+St(i+24*al-1)+St(i+25*al-1)+St(i+26*al-1)+St(i+27*al-1)+St(i+28*al-1)+St(i+29*al-1)); %V0
        dSt(i+31*al,1)=dSt(i+31*al,1) + (1-sc2*v2(round(t),i))*(1-sc3*v3(round(t),i))*u(i-1)*St(i+31*al-1) + sc1*v1(round(t),i)*u(i-1)*(St(i+3*al-1)+St(i+4*al-1)) + (sc2*v2(round(t),i)+sc3*v3(round(t),i))*u(i-1)*(St(i+17*al-1)+St(i+18*al-1)+St(i+30*al-1)); %V1
        dSt(i+32*al,1)=dSt(i+32*al,1) + u(i-1)*St(i+32*al-1) + sc1*v1(round(t),i)*u(i-1)*(St(i+6*al-1)+St(i+7*al-1)+St(i+9*al-1)) + (sc2*v2(round(t),i)+sc3*v3(round(t),i))*u(i-1)*(St(i+20*al-1)+St(i+21*al-1)+St(i+23*al-1)+St(i+31*al-1)); %V2
        
    end
end
