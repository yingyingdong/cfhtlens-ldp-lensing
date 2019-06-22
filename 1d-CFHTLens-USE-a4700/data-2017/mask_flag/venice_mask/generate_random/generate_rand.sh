#num = [343549,64522,193811,87112]
num=871120
field=4
zm=0.45
magl=-21
mask=../1reverse_reg2/output_w${field}.1arcsec.fits
venice -m ${mask} -r -f outside -npart $num -o w${field}random_xy.fits
python generate_rand.py $field $zm $magl
