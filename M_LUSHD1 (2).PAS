
{Апроксимацiя чебишовським наближенням} 
uses crt; 
	var 
	{промiжнi матрицi}     u2:array[1..2,1..2] of double;     u3:array[1..3,1..3] of double; 
	u4:array[1..4,1..4] of double; 
	u5:array[1..5,1..5] of double; 
	function det3:double; {визначник 3-мiрної матрицi u3}     
begin       det3:=u3[1,1]*u3[2,2]*u3[3,3]+u3[3,1]*u3[1,2]*u3[2,3] 
		+u3[2,1]*u3[3,2]*u3[1,3]-u3[3,1]*u3[2,2]*u3[1,3] 
		-u3[2,1]*u3[1,2]*u3[3,3]-u3[1,1]*u3[3,2]*u3[2,3];     
end; 
	function det4:double; {визначник 4-мiрної матрицi u4}     
var k,l:word; aa:double;     
begin 
	  aa:=0; 
	  for k:=1 to 3 do for l:=1 to 3 do u3[k,l]:=u4[k,l]; 
	  aa:=aa+u4[4,4]*det3; 
	  for k:=1 to 3 do u3[k,3]:=u4[k,4]; 
	  aa:=aa-u4[4,3]*det3; 
	  for k:=1 to 3 do u3[k,2]:=u4[k,3]; 
	  aa:=aa+u4[4,2]*det3; 
	  for k:=1 to 3 do u3[k,1]:=u4[k,2]; 
	  aa:=aa-u4[4,1]*det3; 
	  det4:=aa; 
	end; 
	function det5:double; {визначник 5-мiрної матрицi u5}     
var k,l:word; aa:double;     
begin 
	  aa:=0; 
	  for k:=1 to 4 do for l:=1 to 4 do u4[k,l]:=u5[k,l]; 
	  aa:=aa+u5[5,5]*det4; 
	  for k:=1 to 4 do u4[k,4]:=u5[k,5]; 
	  aa:=aa-u5[5,4]*det4; 
	  for k:=1 to 4 do u4[k,3]:=u5[k,4]; 
	  aa:=aa+u5[5,3]*det4; 
	  for k:=1 to 4 do u4[k,2]:=u5[k,3]; 
	  aa:=aa-u5[5,2]*det4; 
	  for k:=1 to 4 do u4[k,1]:=u5[k,2]; 
	  aa:=aa+u5[5,1]*det4; 
	  det5:=aa; 
	end; 
	function sign(x:double):integer; {sign}     
begin 
		if x>0 then sign:=1 else if x<0 then sign:=-1 else sign:=0;     
end; 

{----------- тут завдавати апроксимовану функцiю -------------- } 
	function fun(x:double):double; {функцiя, яка апроксимується}     
begin 
{		fun:=sin(x); {} 
{}		fun:= log(x); {}     
end; 

	function w(x:double):double; {апроксимацiйна вага w}     
begin 
		w:=1; {найкраще абсолютне наближення} 
{		w:=fun(x); {найкраще вiдносне наближення}     
end; 

	var a,b:double; {промiжок} 
	eps, {задана точнiсть} 
	ro,r1, 
	d0,d1,xx, 
	dx, {крок для визначення максимального вiдхилення на промiжку} 
	xm, {точка, в якiй максимальне вiдхилення} 
	d,tt:double; 
	t0,t1:array[0..4] of double; {точки альтернансу: попередня i наступна iтерацiї} 
	nm, {кiлькiсть точок на пормiжку для визначення максимального вiдхилення} 
	i,k,l,m,m1,v:word; {m - порядок полiнома; m=3} 
	aa:array[0..4] of double; { aa[0]-aa[3] - коефiцiєнти шуканого апроксимуючого полiнома, aa[4] - miu} 
	cc:array[1..5,0..4] of double; {коефiцiєнти системи рiвнянь для aa} 
	bb:array[1..5] of double; {вiльнi члени цiєї системи рiвнянь} 
	u50:array[1..5,1..5] of double; 
	j:word; {номер iтерацiї} 

	function P(x:double):double; {апроксимуючий полiном P} 
	var ii:word; PP,xx:double;     
begin 
	  PP:=0; xx:=1; 
	  for ii:=0 to m do begin 
	   PP:=PP+aa[ii]*xx; xx:=xx*x; 
	  end; 
	  P:=PP;     
end; 

	begin {основна програма} 
	m:=3; m1:=m+1; 
	a:=0.1; b:=3.14; 
	eps:=0.001; 
	d:=(b-a)/m1; 
	nm:=20; dx:=(b-a)/nm; 
	{початковi значення точок альтернансу} 
	for k:=0 to m1 do t0[k]:=a+d*k; 
	ro:=1.0; aa[4]:=0; {початковi значення для запуску iтерацiйного циклу} 
	{iтерацiї} 
	j:=0; 
	while ro-abs(aa[4])>eps*abs(aa[4]) do begin {iтерацiйний цикл} 
{	  writeln('iтерацiя ',j,' точки альтернансу '); {} 
{	  for k:=0 to m1 do write(t0[k],' '); writeln; if readkey=' ' then; {} 
		{заповнення матрицi cc} 
	  for i:=1 to m+2 do begin 
		tt:=t0[i-1]; 
		cc[i,0]:=1; 
		for k:=1 to m do cc[i,k]:=cc[i,k-1]*tt; 
		cc[i,4]:=w(tt)*(1.0-2.0*((i+1) mod 2)); 
		bb[i]:=fun(tt); 
	  end; 
		{роздрук матрицi cc} 
{}	  writeln('матриця C '); {} 
{}	  for i:=1 to m+2 do begin {} 
{}		for k:=0 to m1 do write(cc[i,k],' '); writeln; if readkey=' ' then; {} 
{}    end; {} 
		{розв'язування системи рiвнянь методом Крамера для m1+1=5} 
	  for i:=1 to 5 do for k:=1 to 5 do u5[i,k]:=cc[i,k-1]; 
	  for i:=1 to 5 do for k:=1 to 5 do u50[i,k]:=u5[i,k]; 
	  d0:=det5; {основний визначник} 
{	  writeln('d0 ', d0); {} 
	  if d0=0 then begin 
		writeln('Система рiвнянь вироджена, iтерацiя ',j); if readkey=' ' then;       end; 
	  for k:=1 to 5 do begin 
		for i:=1 to 5 do for l:=1 to 5 do u5[i,l]:=u50[i,l]; 
		for i:=1 to 5 do u5[i,k]:=bb[i]; 
		d1:=det5; {додатковi визначники} 
{  write('d1',k,' ',d1,' '); {}
		aa[k-1]:=d1/d0; 
	  end; 
{	  writeln; {} 
{}	  write('Iтерацiя ',j,' розв. '); {} 
{}	  for i:=0 to m1 do write(aa[i],' '); writeln; if readkey=' ' then; {} 
{}	  writeln('fun в точках '); {} 
{}	  for k:=0 to m1 do write(fun(t0[k]),' '); writeln; if readkey=' ' then; {} 
{}	  writeln('P в точках '); {} 
{}	  for k:=0 to m1 do write(P(t0[k]),' '); writeln; if readkey=' ' then; {} 
{}	  writeln('ro в точках '); {} 
{}	  for k:=0 to m1 do write(fun(t0[k])-P(t0[k]),' '); writeln; if readkey=' ' then; {}         {визначення ro} 
	  ro:=0; xm:=0; xx:=a; 
	  while (xx<=b) do begin 
		r1:=abs(fun(xx)-P(xx)); 
		if r1>ro then begin ro:=r1; xm:=xx; end; 
		xx:=xx+dx; 
	  end; 
{	  writeln('ro ',ro,' xm ',xm); if readkey=' ' then;   {} 
{	  writeln('miu ',aa[4]); if readkey=' ' then;       {} 
	  if ro-abs(aa[4])>eps*abs(aa[4]) then begin 
			{змiна точок альтернансу для наступної iтерацiї} 
		if xm<t0[0] then begin {випадок 2} 
{		  writeln('випадок 2 '); if readkey=' ' then; {} 
		  if sign(fun(xm)-P(xm))=sign(fun(t0[0])-P(t0[0])) then t0[0]:=xm 
		  else begin 
			for k:=m1 downto 1 do t0[k]:=t0[k-1]; 
			t0[0]:=xm; 
		  end; 
		end else if xm>t0[m1] then begin {випадок 3} 
{		  writeln('випадок 3 '); if readkey=' ' then; {} 
		  if sign(fun(xm)-P(xm))=sign(fun(t0[m1])-P(t0[m1])) then t0[m1]:=xm 
		  else begin 
			for k:=0 to m do t0[k]:=t0[k+1]; 
			t0[m1]:=xm; 
		  end; 
		end else begin {випадок 1} 
{		  writeln('випадок 1 '); if readkey=' ' then; {} 
		  v:=0; while (v<m+1) and (xm>=t0[v+1]) do v:=v+1; 
{		  writeln('v ',v); if readkey=' ' then; {} 
		  if sign(fun(xm)-P(xm))=sign(fun(t0[v])-P(t0[v])) then t0[v]:=xm 
		  else t0[v+1]:=xm;         end; 
		j:=j+1; 
		
	  end; 
	end; 
		{розв'язок} 
	writeln('Кiлькiсть iтерацiй ',j,' miu ',aa[4],' Результат: '); 
	for i:=0 to m do write('a(',i,')=',aa[i],' '); writeln; if readkey=' ' then; 
  end. 

