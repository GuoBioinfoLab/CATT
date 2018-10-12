int hamming( char a[], char b[], int size )
{
    int i,diff = 0;
    for(i=0;i<size;i++)
    {
        if( a[i] !=  b[i] ) 
        {
            diff++;
            if( diff ==  3 ) return 3; 
        }
        
    }
    return diff;
}

int hammingv2( char a[], char b[], int size, int upper)
{
    int i,diff = 0;
    for(i=0;i<size;i++)
    {
        if( a[i] != b[i] )
        {
            diff++;
            if( diff > upper ) return 0;
        }
    }
    return 1;
}
