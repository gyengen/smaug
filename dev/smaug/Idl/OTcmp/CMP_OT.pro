PRO Initialisation
	print
	print,'_________________________________SAC AND SMAUG comparison: Orszag-Tang test'
	print,'Please feel free to contact me if you need any further information or help.'
	print,'Norbert Gyenge (e-mail: n.g.gyenge@sheffield.ac.uk)'
	spawn, 'rm -rf cmp_result/'
	spawn, 'mkdir -p cmp_result'
END

FUNCTION filaname_sort, string_input
	smaug_out=FILE_SEARCH(string_input)
	print, smaug_out
	order=fltarr(n_elements(smaug_out))
	FOR i=0, n_elements(smaug_out)-1 DO BEGIN
	    theNumber = 0L
	    buff1=STRSPLIT(smaug_out[i], '_',/EXTRACT)
	    buff2=STRSPLIT(buff1[1], '.',/EXTRACT)
	    reads, buff2[0], buff3
	    order[i]=buff3
	ENDFOR
	RETURN, smaug_out[SORT(order)]
END

PRO SAC, B_SAC, sac_t, sac_rho, sac_energy
	headline='                                                                               '
	it=long(1)
	ndim=long(1)
	neqpar=long(1)
	nw=long(1)
	varname='                                                                               '
	time=double(1)
	dum=long(1)
	dumd=long(1)
	nn=0
	mu=1.257E-1
	nn_i=0
	readu, 1, headline
	readu, 1, it, time, ndim, neqpar, nw
	sac_t=time
	ndim=abs(ndim)
	nx=lonarr(ndim)
	readu,1,nx
	eqpar=dblarr(neqpar)
	readu,1,eqpar
	readu,1,varname
	n1=nx(0)
	n2=nx(1)
	x=dblarr(n1,n2,ndim)
	IF (nn EQ 0) THEN w=dblarr(n2,n1,nw)
	wi=dblarr(n1,n2)
	readu,1,x
	FOR iw=0,nw-1 DO BEGIN
		readu,1, wi
		w(*,*,iw)=rotate(wi,1)
	ENDFOR
	B=dblarr(n2,n1)
	B_SAC=sqrt(((w[*,*,8]+w[*,*,4])^2.0+(w[*,*,9]+w[*,*,5])^2.0))
	sac_rho = TOTAL(w(*,*,0)+w(*,*,7)) ; rho + rhoB
	sac_energy = TOTAL(w(*,*,3)) ; energy
END

PRO SMAUG, smaug_filename, B_SMAUG, smaug_t, smaug_rho, smaug_energy
	headline='                                                                               '
	it=long(1)
	ndim=long(1)
	neqpar=long(1)
	nw=long(1)
	varname='                                                                               '
	time=double(1)
	dum=long(1)
	dumd=long(1)
	nn=0
	mu=1.257E-1
	nn_i=0
	key=0
	unit=0
	openr, 2, smaug_filename
	readu, 2, headline
	readu, 2, it, time, ndim, neqpar, nw
	smaug_t=time
	ndim=abs(ndim)
	nx=lonarr(ndim)
	readu,2,nx
	eqpar=dblarr(neqpar)
	readu,2,eqpar
	readu,2,varname
	n1=nx(0)
	n2=nx(1)
	x=dblarr(n1,n2,ndim)
	IF (nn EQ 0) THEN w=dblarr(n2,n1,nw)
	wi=dblarr(n1,n2)
	readu,2,x
	FOR iw=0,nw-1 DO BEGIN
		readu,2, wi
		w(*,*,iw)=rotate(wi,1)
	ENDFOR
	B=dblarr(n2,n1)
	B_SMAUG=sqrt(((w[*,*,8]+w[*,*,4])^2.0+(w[*,*,9]+w[*,*,5])^2.0))
	smaug_rho = TOTAL(w(*,*,0)+w(*,*,7)) ; rho+rhoB
	smaug_energy = TOTAL((w(*,*,3))) ; energy
	close, 2
END

PRO VISUAL, B_sac, B_smaug, iter, sac_t, smaug_t, cmp_energy, cmp_rho, x, number_of_step
	charsize=1
	FILENAME1='cmp_result/frame_'+STRTRIM(string(iter),1)+'.ps'
	FILENAME2='cmp_result/frame_'+STRTRIM(string(iter),1)+'.png'
	cgPS_Open, FILENAME1
	cgDisplay, WID=1
	cgLoadCT, 13, /Brewer
	pos = cgLayout([2,2], OXMargin=[7,2], OYMargin=[5,12], XGap=9, YGap=11)

	;1st subpanel - SAC
	p = pos[*,0]
	cgContour, B_sac[*,*], nlevels=255, /FILL ,Position=p, /NOERASE, xtitle='x', ytitle='z', charsize=charsize
	cgColorBar, position=[p[0], p[3]+0.035, p[2], p[3]+0.055], RANGE=[min(B_sac[*,*]),max(B_sac[*,*])], charsize=charsize

	;2nd subpanel - SMAUG
	p = pos[*,1]
	cgContour, B_smaug[*,*], nlevels=255, /FILL, Position=p, /NOERASE, charsize=charsize, xtitle='x', ytitle='z'
	cgColorBar, position=[p[0], p[3]+0.035, p[2], p[3]+0.055], RANGE=[min(B_sac[*,*]),max(B_sac[*,*])], charsize=charsize

	;3rd subpanel - Diff image
	p = pos[*,2]
	cgContour, abs(B_sac[*,*]-B_smaug[*,*]), nlevels=255, /FILL, Position=p, /NOERASE, charsize=charsize, xtitle='x', ytitle='z'
	cgColorBar, position=[p[0], p[3]+0.035, p[2], p[3]+0.055], RANGE=[min(abs(B_sac[*,*]-B_smaug[*,*])),max(abs(B_sac[*,*]-B_smaug[*,*]))], charsize=charsize, Format='(e6.0)'

	;4th subpanel - diff plot
	p = pos[*,3]
	IF MAX(cmp_rho) GT MAX(cmp_energy) THEN maximum=MAX(cmp_rho)
	IF MAX(cmp_rho) LE MAX(cmp_energy) THEN maximum=MAX(cmp_energy)
	cgplot, x, cmp_energy, Position=p, xtitle='i', ytitle='$\delta$', Color='red', charsize=charsize, YRANGE=[0, maximum+maximum*0.2], XRANGE=[0, number_of_step], /NOERASE
	cgplot, x, cmp_rho, Position=p, Color='dodger blue', charsize=charsize, XRANGE=[0, number_of_step], /Overplot, /NOERASE
	AL_Legend, ['E', '$\rho$+$\rho$$\subb$'], LineStyle=[0,0], Color=['red','dodger blue'], bthick=1, charsize=1 , box=0, vertical=0

	cgText, 0.2, 0.95, /Normal, 'OT test', Alignment=0.5, Charsize=1.2
	cgText, 0.4, 0.95, /Normal, 'i = '+STRTRIM(STRING(iter),1), Alignment=0.5, Charsize=1.2
	cgText, 0.6, 0.95, /Normal, 't = '+STRTRIM(STRMID(STRING(FLOAT(sac_t)), 0, 8),1), Alignment=0.5, Charsize=1.2
	cgText, 0.8, 0.95, /Normal, 'dim = '+STRTRIM(STRING(n_elements(B_sac[*,0])),1)+'x'+STRTRIM(STRING(n_elements(B_sac[0,*])),1), Alignment=0.5, Charsize=1.2
	cgPS_Close
	spawn, 'convert -density 100 '+FILENAME1+' -quality 100 -rotate +90 -flatten +repage '+FILENAME2
	spawn, 'rm '+FILENAME1
END

FUNCTION ERROR_CHECK, B_SAC, B_SMAUG, smaug_t, sac_t, iter
	error=0
	IF iter EQ 0 THEN BEGIN
		IF n_elements(B_SAC[0,*]) NE n_elements(B_SMAUG[0,*]) THEN BEGIN
			print, 'Error: Different dimensions.' 
			help, B_SAC
			help, B_SMAUG
			error=1
		ENDIF
	ENDIF
	IF (smaug_t-sac_t) NE 0 THEN BEGIN
		print, 'Error: not equal time steps.'
		error=1
	ENDIF
	RETURN, error
END

;************************ INPUT FOLDERS *************************
SMAUG_FOLDER='/data/cs1ngg/smaug/out/*.out'
SAC_FOLDER='/data/cs1ngg/sac_working/results/zeroOT.out'

;***************************** MAIN *****************************
Initialisation
smaug_out=filaname_sort(SMAUG_FOLDER) ; Read SMAUG result files and sort them
number_of_step=n_elements(smaug_out) ; number of simulation time steps
openr, 1, SAC_FOLDER, /f77_unf ; open SAC output file
SAC, B_SAC, sac_t ;ignore the first step
iter=0  ; iteration step
cmp_energy=[] ; energy comparison between sac and smaug -2D plot, 4th subpanel
cmp_rho=[] ; density comparison between sac and smaug, 4th subpanel
x=[] ; iteration step array, 4th subpanel
WHILE NOT(eof(1)) DO BEGIN
	SAC, B_SAC, sac_t, sac_rho, sac_energy ;read one layer of SAC
	SMAUG, smaug_out[iter], B_SMAUG, smaug_t, smaug_rho, smaug_energy ;read one layer of SMAUG
	signal=ERROR_CHECK(B_SAC, B_SMAUG, smaug_t, sac_t, iter) ;check dimensions and time 
	IF signal EQ 0 THEN BEGIN 
		x=[x,iter] ;fill up the 2D plot
		cmp_energy=[cmp_energy, ABS(sac_energy-smaug_energy)] ;fill up the 2D plot
		cmp_rho=[cmp_rho, ABS(sac_rho-smaug_rho)] ;fill up the 2D plot
		VISUAL, B_SAC, B_SMAUG, iter, sac_t, smaug_t, cmp_energy, cmp_rho, x, number_of_step ;create the figure
		iter=iter+1
	ENDIF ELSE BEGIN
		print, 'Program aborted.'
		BREAK
	ENDELSE
ENDWHILE
close, 1
END
