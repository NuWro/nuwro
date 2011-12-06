while(<>)
{
	$i++;
	if($i)
	{
	    trim;
	@a=split(/\s+/);
		if($a[2]!=0)
		{
		  print "$a[0] ",$a[1]/$a[2],"\n";
		 }
    }
}
