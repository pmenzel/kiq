#pragma once

#define KIQ_VERSION_MAJOR 0
#define KIQ_VERSION_MINOR 1
#define KIQ_VERSION_PATCH 0

#define KIQ_VERSION_SUFFIX ""

#define KIQ_STR_HELPER(x) #x
#define KIQ_STR(x) KIQ_STR_HELPER(x)

#define KIQ_VERSION_STRING                                                \
	(KIQ_STR(KIQ_VERSION_MAJOR) "." KIQ_STR(KIQ_VERSION_MINOR) "." KIQ_STR( \
	 KIQ_VERSION_PATCH) KIQ_VERSION_SUFFIX)

