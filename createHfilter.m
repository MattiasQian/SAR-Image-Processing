% createHfilter
%
% This is an incomplete script that creates the Hfilter.
%

Ronear=10446;       %Distance to near range [m]
delta_R=4;          %Pixel size in range direction [m]
NOfRows=350;        %Number of rows


lambda=0.0566;      %Wavelength [m]

n=-N/2:N/2;         %Update n after you determine the filter length!


PRF_v=  2.33  ;         %PRF divided by aircraft speed [m^-1]


Hfilter=zeros(NOfRows,length(n));

for row=1:NOfRows
    R= Ronear + (row-1) * delta_R;
    Hfilter(row,:)= exp(1j * 2*pi * n.^2 / (lambda * R * PRF_v^2));
end;

