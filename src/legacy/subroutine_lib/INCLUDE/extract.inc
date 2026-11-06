C Include file for extract_general.f subroutine and related programs



         integer*4 maxofac
	 parameter(maxofac = 16)
	 integer*4 maxchan
	 parameter(maxchan = 4096)
	 integer*4 maxout
	 parameter(maxout = maxchan*maxchan)
	 integer*4 maxnt
	 parameter(maxnt = maxout*maxofac)
