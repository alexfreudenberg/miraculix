XSMALLSIZE=21000
SMALLSIZE=101000
MEDIUMSIZE=750000
LARGESIZE=3100000

SNPSIZE=50241
#SNPSIZE=12803

all:	plink small medium large

#get plink
plink:
	wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip
	unzip plink*zip
	rm LICENSE plink_linux_x86_64_20190617.zip prettify toy.map toy.ped

xsmall:	plink
	head -n $(SNPSIZE) geno.50241  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(XSMALLSIZE) --freq --make-bed --out xsmall
	./generateallelefreq.sh xsmall.frq xsmall.freq
	$(RM) *.frq *.log *.simfreq
	
small:	plink
	head -n $(SNPSIZE) geno.50241  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(SMALLSIZE) --freq --make-bed --out small
	./generateallelefreq.sh small.frq small.freq
	$(RM) *.frq *.log *.simfreq

medium:	plink
	head -n $(SNPSIZE) geno.50241  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(MEDIUMSIZE) --freq --make-bed --out medium
	./generateallelefreq.sh medium.frq medium.freq
	$(RM) *.frq *.log *.simfreq

large:	plink
	head -n $(SNPSIZE) geno.50241  >geno.sim
	./plink --simulate geno.sim --simulate-ncases $(LARGESIZE) --freq --make-bed --out large
	./generateallelefreq.sh large.frq large.freq
	$(RM) *.frq *.log *.simfreq

clean:	cleanplink cleandata

cleanplink:
	rm -fv plink*.zip toy.* prettify plink LICENSE

cleandata:	
	rm -fv *small* *medium* *large* geno.sim
