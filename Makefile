example:
	bld/build.sh nci cice4 25x29
access-om:
	bld/build.sh nci access-om 360x300
access-cm:
	bld/build.sh nci access-cm 360x300
access-hom:
	bld/build.sh nci access-om 1440x1080
access-hcm:
	bld/build.sh nci access-cm 1440x1080

clean:
	rm -rf build_*
