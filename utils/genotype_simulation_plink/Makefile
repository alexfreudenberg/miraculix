XSMALLSIZE=21000
SMALLSIZE=101000
MEDIUMSIZE=750000
LARGESIZE=3100000

SNPSIZE_XSMALL=12803
SNPSIZE_SMALL=50241
SNPSIZE_MEDIUM=250000
SNPSIZE_LARGE=1000000000

all: plink small medium large few_snps medium_snps many_snps

geno.1m:
	julia create_sim_file.jl

#get plink
plink: geno.1m
	wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip
	unzip plink*zip
	rm LICENSE plink_linux_x86_64_20190617.zip prettify toy.map toy.ped

xsmall:	plink
	head -n $(SNPSIZE_XSMALL) geno.1m  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(XSMALLSIZE) --freq --make-bed --out $@
	./generateallelefreq.sh xsmall.frq xsmall.freq
	$(RM) *.frq *.log *.simfreq
	
small:	plink
	head -n $(SNPSIZE_SMALL) geno.1m  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(SMALLSIZE) --freq --make-bed --out $@
	./generateallelefreq.sh small.frq small.freq
	$(RM) *.frq *.log *.simfreq

medium:	plink
	head -n $(SNPSIZE_SMALL) geno.1m  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(MEDIUMSIZE) --freq --make-bed --out $@
	./generateallelefreq.sh medium.frq medium.freq
	$(RM) *.frq *.log *.simfreq

large:	plink
	head -n $(SNPSIZE_SMALL) geno.1m  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(LARGESIZE) --freq --make-bed --out $@
	./generateallelefreq.sh large.frq large.freq
	$(RM) *.frq *.log *.simfreq

few_snps: plink
	head -n $(SNPSIZE_SMALL) geno.1m  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(XSMALLSIZE) --freq --make-bed --out $@
	./generateallelefreq.sh $@.frq $@.freq
	$(RM) *.frq *.log *.simfreq

medium_snps: plink
	head -n $(SNPSIZE_MEDIUM) geno.1m  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(XSMALLSIZE) --freq --make-bed --out $@
	./generateallelefreq.sh $@.frq $@.freq
	$(RM) *.frq *.log *.simfreq

many_snps: plink
	head -n $(SNPSIZE_LARGE) geno.1m  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(XSMALLSIZE) --freq --make-bed --out $@
	./generateallelefreq.sh $@.frq $@.freq
	$(RM) *.frq *.log *.simfreq

clean:	cleanplink cleandata

cleanplink:
	rm -fv plink*.zip toy.* prettify plink LICENSE

cleandata:	
	rm -fv *small* *medium* *large* geno.sim
