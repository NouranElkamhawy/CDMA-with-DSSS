clc
clear;
close all;

%Requirements
%1)Generate the QPSK sequences (±1 ± j) of N = 4 users, each composed of randomly generated K = 5 symbols
N=4;
K=5;
randSig=zeros(N,K);

%QPSK_Inphase=(((randi([0 1],K,1)))*2-1);                      %In-phase symbol generation
%QPSK_Outphase=(((randi([0 1],K,1)))*2-1);                     %Quadrature symbol generation 
for i=1:1:N
       randSig(i,:) = (((randi([0 1],K,1)))*2-1) + 1i*(((randi([0 1],K,1)))*2-1);
end

QPSK_Inphase_real = real(randSig);      %In-phase symbol generation
QPSK_Quadphase_img = imag(randSig);     %Quadrature symbol generation
  
 %Another Method%
 
%%% 4 users first we will generate random array of 1 and -1 real
%%% and another random stream of bits for the imaginary part
% Real_N1=randi([-1 ,1],1,5);
% Real_N1(~Real_N1)=-1; %%%%%%%% to convert any zero to -1 
% Real_N2=randi([-1 ,1],1,5);
% Real_N2(~Real_N2)=-1;
% Real_N3=randi([-1 ,1],1,5);
% Real_N3(~Real_N3)=-1;
% Real_N4=randi([-1 ,1],1,5);
% Real_N4(~Real_N4)=-1;
% Imag_N1=randi([-1 ,1],1,5);
% Imag_N1(~Imag_N1)=-1; %%%%%%%% to convert any zero to -1 
% Imag_N2=randi([-1 ,1],1,5);
% Imag_N2(~Imag_N2)=-1;
% Imag_N3=randi([-1 ,1],1,5);
% Imag_N3(~Imag_N3)=-1;
% Imag_N4=randi([-1 ,1],1,5);
% Imag_N4(~Imag_N4)=-1;
% for i=1:1:5    %%%%%%%%% to generate the QPSK Sequence
% N1(i)=(Real_N1(i)+j*Imag_N1(i));
% N2(i)=(Real_N2(i)+j*Imag_N2(i));
% N3(i)=(Real_N3(i)+j*Imag_N3(i));
% N4(i)=(Real_N4(i)+j*Imag_N4(i));
% end
%  disp(N1)
%  disp(N2)
%  disp(N3)
%  disp(N4)
%  
 %%%
 
 
%2)Generate the maximal length spreading codes for the N users (i.e. N= 4 users)
Maximal_length=7;   %Length: 7, 15 and 63 enter it manually or using 
%another method
%Maximal_length = input("enter maximal length: ");
% while (Maximal_length ~=7 && Maximal_length ~=15 && Maximal_length ~=63)
%     Maximal_length = input("Renter maximal length: ");
% end
if Maximal_length == 7
   % disp("test 7");     %for debug
    PN_Spreaded_Sequence_in=zeros(N,Maximal_length);
    PN_Spreaded_Sequence_quad=zeros(N,Maximal_length);
    Genrator_Polynomial=[3 2 0];
    for i=1:1:N     %%%%
        h_in =commsrc.pn('GenPoly',Genrator_Polynomial,'Shift',i,'NumBitsOut',Maximal_length);
        h_quad=commsrc.pn('GenPoly',Genrator_Polynomial,'Shift',i+N,'NumBitsOut',Maximal_length);
        %Another method: different initial state for each iteration so that
        %we can get shifted
      
        PN_Spreaded_Sequence_in(i,:)=generate(h_in);
        PN_Spreaded_Sequence_quad(i,:)=generate(h_quad);
    end
elseif Maximal_length == 15
   % disp("test 15");    %for debug
    PN_Spreaded_Sequence_in=zeros(N,Maximal_length);
    PN_Spreaded_Sequence_quad=zeros(N,Maximal_length);
    Genrator_Polynomial=[4 3 0];
    for i=1:1:N
        h_in =commsrc.pn('GenPoly',Genrator_Polynomial,'Shift',i,'NumBitsOut',Maximal_length);
        h_quad=commsrc.pn('GenPoly',Genrator_Polynomial,'Shift',i+N,'NumBitsOut',Maximal_length);
        %Another method: different initial state for each iteration so that
        %we can get shifted
        
        PN_Spreaded_Sequence_in(i,:)=generate(h_in);
        PN_Spreaded_Sequence_quad(i,:)=generate(h_quad);
    end
elseif Maximal_length == 63
   % disp("test 63");    %for debug
    PN_Spreaded_Sequence_in=zeros(N,Maximal_length);
    PN_Spreaded_Sequence_quad=zeros(N,Maximal_length);
    Genrator_Polynomial=[6 5 0];
    for i=1:1:N
        h_in =commsrc.pn('GenPoly',Genrator_Polynomial,'Shift',i,'NumBitsOut',Maximal_length);
        h_quad=commsrc.pn('GenPoly',Genrator_Polynomial,'Shift',i+N,'NumBitsOut',Maximal_length);
        %Another method: different initial state for each iteration so that
        %we can get shifted
        
        PN_Spreaded_Sequence_in(i,:)=generate(h_in);
        PN_Spreaded_Sequence_quad(i,:)=generate(h_quad);
    end
end

PN_Spreaded_Sequence_in(PN_Spreaded_Sequence_in==0)=-1;
PN_Spreaded_Sequence_quad(PN_Spreaded_Sequence_quad==0)=-1;


%3) Spread the signal by multiplying each QPSK symbol with the spreading code
Spreaded_Signal_InPhase_real = zeros(N,K*Maximal_length);
Spreaded_Signal_Quad_img = zeros(N,K*Maximal_length);
for i=1:1:N
    %Note that each of the in-phase and quadrature-phase components will be multiplied by a separate spreading code
    %in-phase component
    Spreaded_Signal_InPhase_real(i,:)=kron(QPSK_Inphase_real(i,:),PN_Spreaded_Sequence_in(i,:));
    %quadrature-phase component
    Spreaded_Signal_Quad_img(i,:)=kron(QPSK_Quadphase_img(i,:),PN_Spreaded_Sequence_quad(i,:));
    
end

%4) Convolute the spreaded signal with the channel impulse response
Channel_Case = 1 ; %change it according to channel manually or
%another method
%Channel_Case = input("enter maximal length: ");
% while (Channel_Case ~=1 && Channel_Case ~=2 && Channel_Case ~=3)
%     Maximal_length = input("Renter maximal length: ");
% end

if Channel_Case == 1
    %disp("test 1");
    h=zeros(N,1);
    h(1:N,:)=1; 
elseif Channel_Case == 2
    %disp("test 2");
    h=zeros(N,4);
    h(1,:)=[1 0.8 0.3 0.05];
    h(2,:)=[1 0.8 0.3 0.05];
    h(3,:)=[1 0.8 0.3 0.05];
    h(4,:)=[1 0.8 0.3 0.05];   
elseif Channel_Case == 3
    %disp("test 3");
    h=zeros(N,4);
    h(1,:)=[1 0.8 0.3 0.05];
    h(2,:)=[0.9 0.7 0.4 0.1];
    h(3,:)=[1 0.2 0.01 0];
    h(4,:)=[1 0.2 0.01 0];
end

%Transmitted_Signal_InPhase_real =zeros(N,K*Maximal_length+2);
%Transmitted_Signal_Quad_img =zeros(N,K*Maximal_length+2);
for i=1:1:N
    Transmitted_Signal_InPhase_real(i,:) = conv(Spreaded_Signal_InPhase_real(i,:),h(i,:));
    Transmitted_Signal_Quad_img(i,:) = conv(Spreaded_Signal_Quad_img(i,:),h(i,:));
end
RecievedSignal_InPhase_real=sum(Transmitted_Signal_InPhase_real);
RecievedSignal_Quad_img=sum(Transmitted_Signal_Quad_img);

%5) At the base station, apply correlator for each of the N users %%%%%
DespreadedMessage_InPhase_real =zeros(N,K);
DespreadedMessage_Quad_img =zeros(N,K);

for i=1:1:N
   for j=1:Maximal_length:(Maximal_length*K)
       
       DespreadedSymbol_InPhase_real = RecievedSignal_InPhase_real(j:j+Maximal_length-1).*PN_Spreaded_Sequence_in(i,:);
       DespreadedMessage_InPhase_real(i,ceil(j/Maximal_length)) =sum(DespreadedSymbol_InPhase_real)/Maximal_length;
       
       DespreadedSymbol_Quad_img = RecievedSignal_Quad_img(j:j+Maximal_length-1).*PN_Spreaded_Sequence_quad(i,:);
       DespreadedMessage_Quad_img(i,ceil(j/Maximal_length)) =sum(DespreadedSymbol_Quad_img)/Maximal_length;
        
   end
end

%6) Finally, apply decision criteria to estimate the transmitted QPSK symbols
DecodedMessage_InPhase_real = DespreadedMessage_InPhase_real;
DecodedMessage_Quad_img = DespreadedMessage_Quad_img;

DecodedMessage_InPhase_real(DecodedMessage_InPhase_real>=0)=1;
DecodedMessage_InPhase_real(DecodedMessage_InPhase_real<0)=-1;

DecodedMessage_Quad_img(DecodedMessage_Quad_img>=0) = 1;
DecodedMessage_Quad_img(DecodedMessage_Quad_img<0) = -1;

%Recieved Random signal
Recieved_RandSig = DecodedMessage_InPhase_real + 1i*DecodedMessage_Quad_img;


%%Deliverables (i.e. the plot)
%Plot the received signal (using stairs()) before and after de-spreading 
%Before de-spreading
figure(1);
subplot(2,1,1);
stairs(RecievedSignal_InPhase_real);
title(['Recieved Message in-phase component before Despreading ',' with PN sequence of Length: ',num2str(Maximal_length),'  case: ',num2str(Channel_Case)]);
subplot(2,1,2);
stairs(RecievedSignal_Quad_img);
title(['Recieved Message quadrature-phase component before Despreading ',' with PN sequence of Length: ',num2str(Maximal_length),'  case: ',num2str(Channel_Case)]);

%After de-spreading
figure(2);
for i=1:N
    subplot(N,1,i)
    stairs(DespreadedMessage_InPhase_real(i,:),'-o');
    ylim ([-3 3]);
    title(['Recieved Message after Despreading Despreaded Message in-phase component of user: ',num2str(i),' with PN sequence of Length: ',num2str(Maximal_length),'  case: ',num2str(Channel_Case)]);
%      subplot(N,2,i+1)
%      stairs(Spreaded_Signal_InPhase_real(i,:));
%      title(['Recieved Message after Despreading Despreaded Message quadrature-phase component of user: ',num2str(i),' with PN sequence of Length: ',num2str(Maximal_length),'  case: ',num2str(Channel_Case)]);

end

figure(3);
for i=1:N
    subplot(N,1,i)
    stairs(DespreadedMessage_Quad_img(i,:),'-o');
    ylim ([-3 3]);
    title(['Recieved Message after Despreading Despreaded Message quadrature-phase component of user: ',num2str(i),' with PN sequence of Length: ',num2str(Maximal_length),'  case: ',num2str(Channel_Case)]);
end

%Then plot estimated transmitted information in the de-spread sequence.

figure(4);
for i=1:N
    subplot(N,1,i)
    stairs(DecodedMessage_InPhase_real(i,:),'-o');
    ylim ([-2 2]);
    title(['Estimated infomation in de-spread sequence in-phase component of user: ',num2str(i),' with PN sequence of Length: ',num2str(Maximal_length),'  case: ',num2str(Channel_Case)]);
end

figure(5);
for i=1:N
    subplot(N,1,i)
    stairs(DecodedMessage_Quad_img(i,:),'-o');
    ylim ([-2 2]);
    title(['Estimated infomation in de-spread sequence quadrature-phase component of user: ',num2str(i),' with PN sequence of Length: ',num2str(Maximal_length),'  case: ',num2str(Channel_Case)]);
    
end


figure(6);
for i=1:N
    subplot(N,1,i)
    stairs(DespreadedMessage_InPhase_real(i,:),'-o');
    ylim ([-3 3]);
    hold on 
    title(['In-phase component of user: ',num2str(i),' with PN sequence of Length: ',num2str(Maximal_length),' case: ',num2str(Channel_Case)]);
%      subplot(N,2,i+1)
      stairs(DecodedMessage_InPhase_real(i,:),'-*');
      hold off
      legend('Recieved Meessage after Despread','Estimated transimited information');
end
figure(7);
for i=1:N
    subplot(N,1,i)
    stairs(DespreadedMessage_Quad_img(i,:),'-o');
    ylim ([-3 3]);
    hold on 
    title(['Quadrature-phase component of user: ',num2str(i),' with PN sequence of Length: ',num2str(Maximal_length),' case: ',num2str(Channel_Case)]);
%      subplot(N,2,i+1)
    stairs(DecodedMessage_Quad_img(i,:),'-*');
    hold off
    legend('Recieved Meessage after Despread','Estimated transimited information');
end