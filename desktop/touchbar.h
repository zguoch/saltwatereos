#ifndef _TOUCHBAR_H
#define _TOUCHBAR_H

#include <QtWidgets>
#import <AppKit/AppKit.h>

// This example shows how to create and populate touch bars for Qt applications.
// Two approaches are demonstrated: creating a global touch bar for the entire
// application via the NSApplication delegate, and creating per-window touch bars
// via the NSWindow delegate. Applications may use either or both of these, for example
// to provide global base touch bar with window specific additions. Refer to the
// NSTouchBar documentation for further details.

// The TouchBarProvider class implements the NSTouchBarDelegate protocol, as
// well as app and window delegate protocols.
@interface TouchBarProvider: NSResponder <NSTouchBarDelegate, NSApplicationDelegate, NSWindowDelegate>

@property (strong) NSCustomTouchBarItem *touchBarItem1;
@property (strong) NSCustomTouchBarItem *touchBarItem2;
@property (strong) NSButton *touchBarButton1;
@property (strong) NSButton *touchBarButton2;

@property (strong) NSObject *qtDelegate;

@end

// Create identifiers for two button items.
static NSTouchBarItemIdentifier Button1Identifier = @"com.myapp.Button1Identifier";
static NSTouchBarItemIdentifier Button2Identifier = @"com.myapp.Button2Identifier";

#endif