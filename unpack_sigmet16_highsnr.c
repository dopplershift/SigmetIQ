/* This routine unpacks Sigmet HighSNR 16-bit packed data */

void unpack_highsnr(unsigned short *inData, float *outData, int len)
{
  int iMan, iExp, i; 
	unsigned short *iptr;
	unsigned short uintData;
	float *optr;

/* The following line performs a byteswap of the 16-bit integer.
/* However, since the data should come in Little Endian format, it should
/* be unnecessary. */
/*	uintData = (uintData << 8) | (uintData >> 8); */

iptr = inData;
optr = outData;

for(i = 0; i < len; i++)
		{
		uintData = *iptr++;
		if(uintData & 0xF000)    /* Exponent > 0  */
				{
				iMan = uintData & 0x7FF; /* Mantissa is first 11 bits */
				iExp = (uintData >> 12) & 0x00F; /* Mantissa is last 4 */
				
				/* Check sign bit (bit 12) and adjust mantissa accordingly */
				if(uintData & 0x0800)
						iMan |= 0xFFFFF000;
				else
						iMan |= 0x00000800;
				
				*optr++ = (float)((float)iMan)*((float)(unsigned short)(1<<iExp))/3.3554432E7;
				}
		else					/* Exponent = 0  */
				*optr++ = (float)((float)(((signed short)uintData) << 20)) / 1.759218603E13;
		}
return;
}
