%function Multidays_mixture_model_Fit(MTSBELDG,MTBl)
%
%Code to find the values of P_random P_serial and P_spatial for which the simulated data fit best the experimental data
%Simply put the experimental data matrices MTSBELDG and MTBl as inputs. The 'simulated_data' file should first be generated using the function Multidays_mixture_model_generate(MTSBELDG)


function Multidays_mixture_model_Fit(MTSBELDG,MTBl)

%load simulation matrices
load('simulated_data'); 

%Generate experimental distributions
[GxD_N,VxD_N,NxD_serial,D_NsegDist] = Distributions_days(MTSBELDG,MTBl);

%computation of error matrices
P = 0:2:100;
A = zeros(length(P),length(P));
B = A;
C = A;
D = A;
nDays = size(GxD_N,2);
mean_sqr_error = zeros(length(P),length(P),10,nDays);
min_error_V = nan(nDays,10);
I_random = zeros(nDays,10);
I_serial = zeros(nDays,10);
for dd = 1:nDays
    for kk = 1:10
        for ii = 1:length(P)
            for jj = 1:length(P)

                G = squeeze(Gx10xPRxPSxD(:,kk,ii,jj,dd));
                A(ii,jj) = mean((G - GxD_N(:,dd)).^2);

                NDV = squeeze(NDVx10xPRxPSxD(:,kk,ii,jj,dd));
                B(ii,jj) = mean((NDV - VxD_N(:,dd)).^2);

                NSB = squeeze(NSBx10xPRxPSxD(:,kk,ii,jj,dd));
                C(ii,jj) = mean((NSB - NxD_serial(:,dd)).^2);

                Nseg = squeeze(Nsegx10xPRxPSxD(:,kk,ii,jj,dd));
                D(ii,jj) = mean((Nseg - D_NsegDist(:,dd)).^2);

            end
        end
        
        [min_V,min_I] = min(A+B+C+D,[],'all');
        [I,J] = ind2sub(size(A),min_I);
        I_random(dd,kk) = I;
        I_serial(dd,kk) = J;
        min_error_V(dd,kk) = min_V;
        mean_sqr_error(:,:,kk,dd) = A+B+C+D;
    end
end

mse_mat = zeros(length(P),length(P),nDays);
Prandom = zeros(nDays,1);
Pserial = zeros(nDays,1);
for dd = 1:nDays
    mse_mat(:,:,dd) = mean(mean_sqr_error(:,:,:,dd),3);
    Prandom(dd) = mean(P(I_random(dd,:)));
    Pserial(dd) = mean(P(I_serial(dd,:)));
end
ndx_random = round(Prandom/2);
ndx_serial = round(Pserial/2);

% Plot m.s.e matrix and overlays of experimental and simulated distributions
figure; 
for dd = 1:nDays
    subplot(5,nDays,dd);    imagesc(mse_mat(:,:,dd).^0.1); axis off
    subplot(5,nDays,dd+nDays); plot(-12:1:12,GxD_N(:,dd),'b',-12:1:12,mean(Gx10xPRxPSxD(:,:,ndx_random(dd),ndx_serial(dd),dd),2),'r');ylim([0 30])
    subplot(5,nDays,dd+2*nDays); plot(-11:1:12,VxD_N(:,dd),'b',-11:1:12,mean(NDVx10xPRxPSxD(:,:,ndx_random(dd),ndx_serial(dd),dd),2),'r');ylim([0 20])
    subplot(5,nDays,dd+3*nDays); plot(1:size(NxD_serial,1),NxD_serial(:,dd),'bo-',1:size(NSBx10xPRxPSxD,1),mean(NSBx10xPRxPSxD(:,:,ndx_random(dd),ndx_serial(dd),dd),2),'r');ylim([0 80]);xlim([0 10])
    subplot(5,nDays,dd+4*nDays); plot(1:size(D_NsegDist,1),D_NsegDist(:,dd),'b',1:size(Nsegx10xPRxPSxD,1),mean(Nsegx10xPRxPSxD(:,:,ndx_random(dd),ndx_serial(dd),dd),2),'r');xlim([0 50])
end

%Plot strategy evolution over time
figure; 
bar([Prandom Pserial 100-(Prandom+Pserial)],1,'stacked','EdgeColor','none')
xlabel('Days');
ylabel('Strategy proportion')

