function result = run_edge(imagefile, dilate, min_rows, min_cols, ...
			    canny_level, filter, direct)

if nargin < 2
  dilate = 5;
end

if nargin < 3
  min_rows=0;
  min_cols=0;
end

if nargin < 5
  canny_level=3;
end

if nargin < 6
  filter=1;
end

if nargin < 7
  direct = '.';
end

I = imread(imagefile);
  
Z=get_sparse_representation(I, 0, canny_level, dilate);


  [r c] = size(Z);

  if r < min_rows
    Z2 = zeros(min_rows,c);

    diff = min_rows - r;

    Z2(ceil(diff/2):ceil(diff/2)+r-1, :) = Z;
    Z=Z2;
    r=min_rows;
  end

  if c < min_cols
    Z2 = zeros(r,min_cols);

    diff = min_cols - c;

    Z2(:, ceil(diff/2):ceil(diff/2)+c-1) = Z;
    Z=Z2;
    c=min_cols;
  end


  Z=uint8(Z);

  file = sprintf('%s/%s.dat', direct, imagefile);
  
  fp=fopen(file, 'wb');
  
  fwrite(fp, size(Z), 'integer*4');

  fwrite(fp, Z, 'uint8');

  fclose(fp);
  
  
