clear all
im_size = 128;
% Square in the middle
square_size = 28;
im_square = zeros(im_size, im_size);
low_ind = (im_size-square_size)/2;
high_ind = (im_size+square_size)/2;
for i = low_ind:high_ind
im_square( low_ind:high_ind, i) = 1;
im_square( i, low_ind:high_ind) = 1;
end
% Load Sheep-Logan image

im_shepplogan = phantom('Modified Shepp-Logan', im_size);



%%theta = input(prompt)           %projection angle input by the user is stored as theta
prompt = 'What is the number of detectors? ';
N_detector= input(prompt) ;          %number of detectors is specified by the user
prompt = 'What is the angle step size? ';
step= input(prompt) ;   
prompt = 'What is the frequency span? ';
A= input(prompt) ;  
          
t_max=(im_size)*sqrt(2);            %largest cross-section of the image
t_range=1:t_max/N_detector:t_max;   %different rays of the same theta are distinguished by their t values
x_range=0:im_size;         %integer values of x         
y_range=0:im_size;   %integer values of y

m=1;
for theta=1:step:180;
k=1;

for t=t_range;
    for x=x_range;                                 
        y_cut(x+1)=(t-x*cosd(theta))/sind(theta);      %y-values corresponding to the integer values of x
        
    end
    points_1=[x_range;y_cut];    
                            %a set of points intersecting the grid
    for y=y_range;
        x_cut(y+1)=(t-y*sind(theta))/cosd(theta);
    end
   points_2=[x_cut;y_range];      %another set of points intersecting the grid
   grid_cuts=[points_1 points_2].';   %complete set of points where one ray intersects our grid+some artifacts
    grid_cuts(isnan(grid_cuts)) = 0;   %artifact1:NaN values which are not useful are taken to be zero
    grid_cuts(grid_cuts<0)=0;        %artifact2:the points which are actually out of the grid are taken to be zero
    grid_cuts(grid_cuts>im_size)=0;  
   grid_cuts=grid_cuts(all(grid_cuts(:,1:2),2),:); %the artifacts described above are eliminated
 
  [Y,I]=sort(grid_cuts(:,1));  %sorting our points according to x
  grid_cuts=grid_cuts(I,:);
  
   image_ind=ceil(grid_cuts);

   
   
   
   line1=0;
   line2=0; %line integrals which will be the value of the projection function at a specific t is initialized
   for i=1:(length(grid_cuts)-2);
   line1=line1+sqrt((grid_cuts(i+1,1)-grid_cuts(i,1))^2+(grid_cuts(i+1,2)-grid_cuts(i,2))^2)*im_square(image_ind(i,1),image_ind(i,2));
   line2=line2+sqrt((grid_cuts(i+1,1)-grid_cuts(i,1))^2+(grid_cuts(i+1,2)-grid_cuts(i,2))^2)*im_shepplogan(image_ind(i,1),image_ind(i,2));
   end
   p_t1(1,k,m)=line1; %results for the im_square
   
   p_t2(1,k,m)=line2; %results for the im_shepplogan
   
   k=k+1;
  
end
m=m+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%RECONSTRUCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_square1 = zeros(im_size, im_size);
im_shepplogan1 = zeros(im_size, im_size);
im_square2 = zeros(im_size, im_size);
im_shepplogan2 = zeros(im_size, im_size);
m=1;

for theta=1:step:180;
  k=1;

fftp=fft(p_t1(1,:,m));
 freqs=0:length(fftp)/(length(fftp)-1):length(fftp); %ramp window

                                     
    p=abs(ifft(fftp.*(A*freqs)));    %filtered backprojections
    
 
fftp2=fft(p_t2(1,:,m));
freqs=0:length(fftp2)/(length(fftp2)-1):length(fftp2);

  p2=abs(ifft(fftp2.*(A*freqs)));
    
 
for t=t_range;
    for x=x_range;                                 
        y_cut(x+1)=(t-x*cosd(theta))/sind(theta);      %y-values corresponding to the integer values of x
        
    end
    points_1=[x_range;y_cut];    
                            %a set of points intersecting the grid
    for y=y_range;
        x_cut(y+1)=(t-y*sind(theta))/cosd(theta);
    end
   points_2=[x_cut;y_range];      %another set of points intersecting the grid
   grid_cuts=[points_1 points_2].';   %complete set of points where one ray intersects our grid+some artifacts
    grid_cuts(isnan(grid_cuts)) = 0;   %artifact1:NaN values which are not useful are taken to be zero
    grid_cuts(grid_cuts<0)=0;        %artifact2:the points which are actually out of the grid are taken to be zero
    grid_cuts(grid_cuts>im_size)=0;  
   grid_cuts=grid_cuts(all(grid_cuts(:,1:2),2),:); %the artifacts described above are eliminated
 
  [Y,I]=sort(grid_cuts(:,1));  %sorting our points according to x
  grid_cuts=grid_cuts(I,:);
  
   image_ind=ceil(grid_cuts);
    
  
   
    for i=1:(length(image_ind)-1);     
   im_square1(image_ind(i,1),image_ind(i,2))= im_square1(image_ind(i,1),image_ind(i,2))+p_t1(1,k,m);
    im_shepplogan1(image_ind(i,1),image_ind(i,2))= im_shepplogan1(image_ind(i,1),image_ind(i,2))+p_t2(1,k,m);

     end
     for i=1:(length(image_ind)-1);     
   im_square2(image_ind(i,1),image_ind(i,2))= im_square2(image_ind(i,1),image_ind(i,2))+p(k);
    im_shepplogan2(image_ind(i,1),image_ind(i,2))= im_shepplogan2(image_ind(i,1),image_ind(i,2))+p2(k);

     end
     
     k=k+1;
     
end

m=m+1;
    end




subplot(2,2,1); imagesc(im_square1); colormap gray; colorbar
title({['square image simple backprojection']; [ 'with ' num2str(N_detector) ' detectors and angle step size '  num2str(step) ]});

subplot(2,2,2); imagesc(im_shepplogan1); colormap gray; colorbar

title({['shepplogan simple backprojection']; [ 'with ' num2str(N_detector) ' detectors and angle step size '  num2str(step)]});

subplot(2,2,3); imagesc(im_square2); colormap gray; colorbar
title({['square image filtered backprojection']; [ 'with ' num2str(N_detector) ' detectors and angle step size '  num2str(step) ]});

subplot(2,2,4); imagesc(im_shepplogan2); colormap gray; colorbar

title({['shepplogan filtered backprojection']; [ 'with ' num2str(N_detector) ' detectors and angle step size '  num2str(step)]});








    

