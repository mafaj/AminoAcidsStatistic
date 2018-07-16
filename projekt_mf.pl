#!/usr/bin/perl
use Math::Complex;
use strict;
#program 2
my $plik_fasta = $ARGV[0];
my $plik_temp="statystyki.csv";
my @array;
my @arr;
my @rozmiar;
my @Amino_cap= qw(- A R N D C Q E G H I L K M F P S T W Y V);
my @amino_small = qw(- a r n d c q e g h i l k m f p s t w y v);
my @amino_names = qw(gap_- alanina_A arginina_R asparagina_N kw_asparginowy_D cysteina_C glutamina_Q kw_glutaminowy_E glicyna_G histydyna_H izoleucyna_I leucyna_L lizyna_K metionina_M fenyloalanina_F prolina_P seryna_S treonina_T tryptofan_W tyrozyna_Y walinia_V);
my @sign;
my $tab = scalar @Amino_cap;
my $nr_seq=0;
if($ARGV[0])
{
	open(FH, $plik_fasta);
	print "Plik $plik_fasta jest otwarty do odczytu\n";
	open(FK,">$plik_temp");
	print "Plik $plik_temp jest otwarty do zapisu\n\nUWAGA do odczytania pliku $plik_temp nalezy uzyc arkusza kalkulacyjnego z ustawieniem separatora na przecinek.\n";
	

	my @a;
	my @name; 
#z pliku fasta przepisuje sekwencje do tablicy i zapisuje w pliku .csv
	while(<FH>)
	{ 	
		if($_ =~ />/)
		{ 
			chomp $_;
			push @name, $_;#zapisuje nazwę sekwencji
		}
		else
		{ 
			push @array, $_;#tablica array przechowuje sekwencje
			
			@a = split('', $_);
			my $d = scalar @a;
			push @rozmiar, ($d);
			@a = join (',',@a); 
			push @arr, @a;  
			$nr_seq = $nr_seq+1;
		}
	}
##	#odnajduję rozmiar sekwencji
	my @size = sort {$b <=> $a} @rozmiar;
	my $max_size=$size[0]-1;
	my $v;
	for(-1..$max_size)
	{
		my $z= join (',',  $_);
		print FK "$z,";#zapisuje do pliku csv nr. kolumn (pozycji)
	}
	print FK "\n";
	
	for(my $m =0; $m<$nr_seq; ++$m)#zapisuje do pliku csv nr.,nazwę i sekwencje
	{
		my $n=$m+1;
		print FK "$n, $name[$m], $arr[$m]";
	}	
	

	my $kolumna = 0;


#ile razy dany aminokwas lub gap w danej kolumnie się pojawił
	for(my $q=0; $q < $tab; ++$q)
	{
		my @gaps;
		foreach(@array)
		{ 
			@a = split('', $_);
			for( $kolumna=0;$kolumna<($max_size); ++$kolumna)
				{ 
					if($a[$kolumna] =~ $Amino_cap[$q] || $a[$kolumna] =~  $amino_small[$q])
					{
						$gaps[$kolumna]=$gaps[$kolumna]+1;
					}
				}
		}
		if($q==0) {@sign = @gaps;}
		@gaps = join(',',@gaps);
		print FK "\n,$amino_names[$q],@gaps";
	}
	
#ile w każdej kolumnie jest aminokwasów w ogóle, wliczając '-'
	print FK "\n\n,--ilosc-z gap'ami [%]--, ";
	for(my $q=0; $q < $tab; ++$q)
	{ 
		print FK "\n,$amino_names[$q],";
		
		for( $kolumna=0;$kolumna<($max_size); ++$kolumna)
		{ 
			my $gap=0;
			my $ilosc_amino=0;
			foreach(@array)
			{
				 @a = split('', $_);
				if($a[$kolumna] =~ $Amino_cap[$q] || $a[$kolumna] =~  $amino_small[$q])
				{
					++$gap;
				}
			}
				my $pro= $gap/$nr_seq*100;
				$pro = eval sprintf('%.2f', "$pro");
				print FK $pro . ",";
		}
	}

#ile w każdej kolumnie jest aminokwasów w ogóle, nie licząc '-'
	print FK "\n\n,-- ilosc-bez gap'ow [%]--, ";
	for(my $q=1; $q < $tab; ++$q)
	{ 
		print FK "\n,$amino_names[$q],"; 
		for( $kolumna=0;$kolumna<($max_size); ++$kolumna)
		{
			my $gap=0;
			my $a=0;
			my $ilosc_amino=0;
			
			foreach(@array)
			{
				 @a = split('', $_);
				if($a[$kolumna] =~ $Amino_cap[$q] || $a[$kolumna] =~  $amino_small[$q])
				{
					++$gap;
				}
			}
				
			my $pro2 =$gap/($nr_seq-$sign[$kolumna])*100;
			my $num = eval sprintf('%.2f', "$pro2");
			print FK $num . ",";
		}
	}
	
	#ile w każdej kolumnie jest różnych aminokwasów (
	#na podstawie identyczności w kolumnach liczę średnią i odchylenie standardowe
	print FK "\n\n,ilosc roznych aminokwasów,";

	my @ilosc_bg;
	my $srednia_identycznosc;
	my @srednia_tab;#wartości do zapisaniaw pliku wynikowym
	my @od_stan_tab; #wartości do zapisaniaw pliku wynikowym
	for($kolumna=0;$kolumna<$max_size;++$kolumna)
	{
		my $w;
		my $gap;
		my @nr_a=();
		my $nr=0;
		my @kolumna_a; #zawiera aminokwasy z jednej kolumny
		for(my $i=0;$i<$nr_seq;++$i)
		{ 
			@a=split('', @array[$i]);
			push @nr_a, $a[$kolumna];
			push @kolumna_a, $a[$kolumna];
			
			if($a[$kolumna] =~ /-/){$gap=1;}
			$nr =  @nr_a;
			for(my $j=0;$j<$nr;++$j)
			{
				my $o=$nr-1;
				if($nr_a[($o)]eq$nr_a[$j] && $o!=$j )
				{
					pop (@nr_a); 
					$nr =  @nr_a;
				}
			}
		}

		#obliczanie identycznosci
		
		my @identyczne;#wyniki poszczegulnych porównań parami
		my $ilosc_porownan;
		my $srednia;
		my $wariancja;
		my $od_stan;
		for(my $k=0;$k<@kolumna_a;++$k)
		{
			for(my $l=0;$l<@kolumna_a;++$l)
			{
				my $el;
				if($k < $l)
				{
					#porównuję każdą parę aminokwasów w kolumnie, tylko raz
					#przyjmuję że para '-' i '-' nie jest identyczna
					
					if($kolumna_a[$k] ne $amino_small[0] && ($kolumna_a[$k] =~ $kolumna_a[$l]))
					{
						$el=1;
						push @identyczne, $el;
						++$srednia;
					}	
					else 
					{ 
						$el=0;
						push @identyczne,$el;
					}
				}
			}
		}
		$ilosc_porownan=@identyczne;
		$srednia/=$ilosc_porownan;
		$srednia_identycznosc+=$srednia;
		
		foreach(@identyczne)
		{
			$wariancja+=($_-$srednia)*($_-$srednia)
		}
		
		$wariancja/=$ilosc_porownan;
		$od_stan = sqrt($wariancja);
		
		$srednia_tab[$kolumna]=eval sprintf('%.2f', "$srednia");
		$od_stan_tab[$kolumna]=eval sprintf('%.2f', "$od_stan");
		
		$ilosc_bg[$kolumna]=$nr-$gap;
		print FK "$nr,";
		
	}
		
		
		print FK "\n,ilosc bez gap'ow,", join(',',@ilosc_bg);
		
		@srednia_tab = join(',',@srednia_tab);
		print FK "\n\n,srednia,@srednia_tab";
		
		print FK "\n,odchylenie standardowe,", join(',',@od_stan_tab);
		print "\nUWAGA program przyjmuje, ze para '-' i '-' nie jest identyczna\n";
		$srednia_identycznosc/=$max_size;
		$srednia_identycznosc/=$max_size;
		$srednia_identycznosc*=100;
		print "\nSrednia identycznosc wszystkich kolumn wynosi: $srednia_identycznosc%";
#------------------------------------*****************************************----------------------------------------
#statystyki dla grupy zdefiniowanej przez użytkownika
	my $odp;
	print "\nCzy chcesz policzyc statystyki dla wlasnej zdefiniowanej grupy? 't'/'n'\n";
	$odp=<STDIN>;
	if($odp =~ /n/|| $odp =~ /N/){print "to nie";}
	elsif ($odp =~ /t/ || $odp =~ /T/)
	{
		
		print "Podaj grupe aminokwasow (WIELKIMI LITERAMI): \n";
		my @grupa;
		my $mino=<STDIN>;
		chomp;
		@grupa= split ('',$mino);
		pop (@grupa);
		my $y;
		my $i;
		my $z = @grupa."";
		
		foreach($y=0;$y<$z;++$y)
		{
			my $licz=0;
			for($i=0;$i<$tab;$i++)
			{
				if($grupa[$y]=~$Amino_cap[$i])
				{
					$licz=1;
				}
			}
			if($licz==/0/){print "taki aminokwas \"$grupa[$y]\" nieistnieje\npopraw grupe\n";}
		}
		
		
		print FK "\n\ngrupa:,", @grupa;
		print FK "\n, % pojawienia sie grupy,"; 
		
		my @score;
		my @procent_grupa=();
		my @obecnosc;
		my $seq;
		#dla każdej kolumny 
			for($kolumna=0;$kolumna<$max_size;++$kolumna)
			{	
				my $gr=0;
				#dla każdej sekwencji 
				
				for($i=0;$i<$nr_seq;++$i)
				{push @obecnosc, 0;
					
					my $kolumna1=$kolumna;
					my $gr1;
					@a=split('', @array[$i]);
					
					#sprawdzam czy pierwszy aminokwas z grupy(podanej przez użyt.) == odczytanemu z sekwencji aminokwasowi
					if($a[$kolumna] =~ $grupa[0])
					{
						
						#jeśli pasuje, sprawdzam czy drugi aminokwas z grupy znajduje się w tej samej sekwencji w kolumnie obok
						++$kolumna1;
						++$gr1;
						
						for($gr1..$z)
						{
							#i dalej czy następne aminokwasy z grupy w następnych kolumnach ($z- ilość amninokwasów w grupie)
							if($a[$kolumna1] =~ $grupa[$gr1])
							{
								++$kolumna1;
								++$gr1;
							}
								#jeśli iterator aminokwasów w grupie(gr1) == ilość aminokwasów w grupie => cała grupa znajduje się w sekwencji
								if($gr1 eq $z)
								{
									for(my $u=$kolumna;$u<$kolumna1;++$u)
									{
										#w odpowiednich kolumnach zwiększ @score
										++$score[$u];
										++$seq;
										++$obecnosc[$i];
									}
								}
						}
					}#jeśli w następnych kolumnach nie ma aminokwasów opowiadajęcym tym w grupie, wracamy do miejsca, w któym znaleźliśmy pierwszwy pasujący aminokwas
				}
				$procent_grupa[$kolumna] = $score[$kolumna]/$nr_seq*100;
				my $p = eval sprintf('%.1f', "$procent_grupa[$kolumna]");
				print FK $p . ",";
			}
			print FK "\n,ilosc pojawienia sie motywu,",join(',', @score);
			print "\n szukana sekwencja '@grupa' wystapila $seq razy w zestawieniu";
			
			my $licznik=0;#liczy ilość sekwencji, w których występuje grupa
			for(my $aa=0;$aa<$nr_seq;++$aa)
			{
				if($obecnosc[$aa]!= /0/)
				{
					++$licznik;
				}
			}
			print "\nNa $nr_seq sekwencji szukany motyw pojawil sie w $licznik, to jest ", $licznik/$nr_seq*100,"%;";
			
			
			print "\nwyswietl nazwy sekwencji w ktorych motyw: \n[1]wystapil \n[2]nie wystapil \n[0]nic nie wyswietlaj\n";
			$odp=<STDIN>;
			if($odp =~ /1/)
			{
				for(my $aa=0;$aa<$nr_seq;++$aa)
				{
					if($obecnosc[$aa] != /0/)
					{
						print $name[$aa],"\n";
					}
				}
			}
			elsif($odp =~ /2/)
			{
				for(my $aa=0;$aa<$nr_seq;++$aa)
				{
					if($obecnosc[$aa] =~ /0/)
					{
						print $name[$aa],"\n";
					}
				}
			}
	}
	else {print "nie wlasciwe polecenie\n";}

}
else
{
	print "Need to get fasta file on the command line\n";
}
close(FH);
close(FK);
exit;