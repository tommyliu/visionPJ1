function result=get_sparse_representation(I, edge_thresh, ...
						     new_stuff, dilate_pixels)
global CACHE_FILE_PREFIX

if nargin < 3
  new_stuff = 0;
end



I_no = I;
%I = get_cached_image(I);

%[r c] = size(I);

%result = zeros(r,c,4);

if new_stuff==1 || new_stuff==0

  f{1} = [-1; 0; 1];
  f{2} = [0 0 1; 0 0 0; -1 0 0];
  f{3} = [-1 0 1];
  f{4} = [-1 0 0; 0 0 0; 0 0 1];

  for i=1:4
    %  result(:,:,i) = abs(conv2(I, f{i}, 'same')) > edge_thresh;
    result(:,:,i) = abs(conv2(I, f{i}, 'same')); %> edge_thresh;
  end

  result=result/255;

  norm = sqrt(sum(result .* result, 3));


  for i=1:4
    result(:,:,i) = result(:,:,i) ./ (norm + edge_thresh); 
  end

  if(new_stuff == 1)
    mask=imdilate(edge(I,'canny'),ones(2,2));

    for i=1:4
      result(:,:,i) = result(:,:,i) .* mask;
    end
  end

elseif new_stuff==2

  H=filter2([-1 0 1], I, 'same');
  V=filter2([-1;0;1], I, 'same');
  warning off MATLAB:divideByZero;
  D=atan(H./V);
  warning on MATLAB:divideByZero;
  D(find(V==0)) = pi/2;

  [aa bb] = size(H);
  Q = zeros(aa,bb);

  Q(find((D < -3*pi/8))) = 1;
  Q(find((D > 3*pi/8))) = 1;
  Q(find((D > -3*pi/8) .* (D < -pi/8))) = 2;
  Q(find((D > -pi/8) .* (D < pi/8))) = 3;
  Q(find((D > pi/8) .* (D < 3*pi/8))) = 4;
  E = edge(I,'canny');
  result = Q .* E;

elseif new_stuff == 3

  if(prod(size(I_no))) == 1
    fname = sprintf('%s/edges/%d.dat',CACHE_FILE_PREFIX , I_no);
    fid=fopen(fname, 'rb');
fname
    if fid ~= -1
       r = fread(fid, [1 1], 'integer*4');
       c = fread(fid, [1 1], 'integer*4');
       result = fread(fid, [r c], 'uint8');
       fclose(fid);
 
       return;
    end

    I=get_cached_image(I);

  end


  H=filter2([-1 0 1], I, 'same');
  V=filter2([-1;0;1], I, 'same');
  warning off MATLAB:divideByZero;
  D=atan(H./V);
  warning on MATLAB:divideByZero;
  D(find(V==0)) = pi/2;

  [aa bb] = size(H);
  Q = zeros(aa,bb);

  Q(find((D < -3*pi/8))) = 1;
  Q(find((D > 3*pi/8))) = 1;
  Q(find((D > -3*pi/8) .* (D < -pi/8))) = 2;
  Q(find((D > -pi/8) .* (D < pi/8))) = 3;
  Q(find((D > pi/8) .* (D < 3*pi/8))) = 4;

  E = edge(I,'canny');
  
  result = Q .* E;

  result2 = zeros(size(result));

  for i=1:4
    result2(:,:,i) = (result == i);

    result2(:,:,i) = imdilate(result2(:,:,i), ones(dilate_pixels));
  end

  result = result2(:,:,1) + result2(:,:,2)*2 + result2(:,:,3)*4 + ...
	   result2(:,:,4)*8;

elseif new_stuff == 4

  if(prod(size(I))==1)
    I=get_cached_image(I);
  end

  H=filter2([-1 0 1], I, 'same');
  V=filter2([-1;0;1], I, 'same');
  warning off MATLAB:divideByZero;
  D=atan(H./V);
  warning on MATLAB:divideByZero;
  D(find(V==0)) = pi/2;

  [aa bb] = size(H);
  Q = zeros(aa,bb);

  Q(find((D < -3*pi/8))) = 1;
  Q(find((D > 3*pi/8))) = 1;
  Q(find((D > -3*pi/8) .* (D < -pi/8))) = 2;
  Q(find((D > -pi/8) .* (D < pi/8))) = 3;
  Q(find((D > pi/8) .* (D < 3*pi/8))) = 4;

  E = edge(I,'canny');
  
  result = Q .* E;

  result2 = zeros(size(result));

  for i=1:4
    result2(:,:,i) = (result == i);
  end

  result3(:,:,1) = result2(:,:,1) | result2(:,:,2) | result2(:,:,4);
  result3(:,:,2) = result2(:,:,1) | result2(:,:,2) | result2(:,:,3);
  result3(:,:,3) = result2(:,:,2) | result2(:,:,3) | result2(:,:,4);
  result3(:,:,4) = result2(:,:,1) | result2(:,:,3) | result2(:,:,4);
  
  result2=result3;

  for i=1:4
    result2(:,:,i) = imdilate(result2(:,:,i), strel('disk',ceil(dilate_pixels/2),0));
  end

  result = result2(:,:,1) + result2(:,:,2)*2 + result2(:,:,3)*4 + ...
	   result2(:,:,4)*8;


end



