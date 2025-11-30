
clear; close all; clc; 

paraPatchSize = 23;

GaborFilters =  getGabor(paraPatchSize,paraPatchSize);

OrthGF = cell(1,5);
for ii=1:5
    OrthGF{ii}=orth(GaborFilters((ii-1)*8+1:ii*8,:)');
end

figure;
subplot(1,2,1);
display_network(GaborFilters',[],[],8);   

subplot(1,2,2);
display_network([OrthGF{1},OrthGF{2},OrthGF{3},OrthGF{4},OrthGF{5}],[],[],8);    