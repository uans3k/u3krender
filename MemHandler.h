#pragma once
#include "CustomDefine.h"

inline void memSet16(void *dest, uint16 data, int count)
{
	// this function fills or sets unsigned 16-bit aligned memory
	// count is number of words

	_asm
	{
		mov edi, dest; //edi points to destination memory
		mov ecx, count; //number of 16 - bit words to move
		mov ax, data; //16 - bit data
		rep stosw; //move data
	}

}

inline void memSet32(void *dest, uint32 data, int count)
{
	// this function fills or sets unsigned 32-bit aligned memory
	// count is number of quads

	_asm
	{
		mov edi, dest; //edi points to destination memory
		mov ecx, count; //number of 32 - bit words to move
		mov eax, data;// 32 - bit data
		rep stosd; //move data
	}

}