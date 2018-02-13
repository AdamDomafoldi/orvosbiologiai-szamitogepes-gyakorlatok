% simple vector
v = [1 2 3 4 5];

% two series
s1 = 1:9;
s2 = 5:10;

% transpose
w = (1:9)';

% column vector
cv = [1;2;3;4];

% lenght of the vector
length(cv);
length(s1);

% random permutation -p1 = length
randperm(4);

% easy arithmetics
a = [1 2 3];
b = [9 8 7]; 
2*a;
a+b;

% index matrix and change the indexed value
a(1:2)=5;
a(1:2)=[4 4];

% basic functions
sum(a);
mean(a);
min(a);
max(a);

% search -> result (if there is any) is the index of the searched value in the matrix
find(b==8);

% string methods  
a='hello';

%concatting strings
strcat('Adam','Donkey');

% cast numeric to string
num2str(pi);

% explicit write out aside the lack of semicolon at the end of the line
name = 'Alice';   
age = 12;
X = [name,' will be ',num2str(age),' this year.'];
disp(X);

% plot
%X = 0:0.1:20; %from zero to 20 with the rate of 0.1 
%figure; hold on; %hold the windows for data
%plot(X, 1+cos(X));
%plot(X, 1+cos(2*X), '-.r'); % set red, dashed, dot line
%legend('1+cos(x)', '1+cos(2x)'); % legend for the curves


 







