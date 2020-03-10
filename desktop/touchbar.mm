#include "touchbar.h"

@implementation TouchBarProvider

- (NSTouchBar *)makeTouchBar
{
    // Create the touch bar with this instance as its delegate
    NSTouchBar *bar = [[NSTouchBar alloc] init];
    bar.delegate = self;

    // Add touch bar items: first, the very important emoji picker, followed
    // by two buttons. Note that no further handling of the emoji picker
    // is needed (emojii are automatically routed to any active text edit). Button
    // actions handlers are set up in makeItemForIdentifier below.
    bar.defaultItemIdentifiers = @[NSTouchBarItemIdentifierCharacterPicker,
                                   Button1Identifier, Button2Identifier];

    return bar;
}

- (NSTouchBarItem *)touchBar:(NSTouchBar *)touchBar makeItemForIdentifier:(NSTouchBarItemIdentifier)identifier
{
    Q_UNUSED(touchBar);

    // Create touch bar items as NSCustomTouchBarItems which can contain any NSView.
    if ([identifier isEqualToString:Button1Identifier]) {
        QString title = "B1";
        self.touchBarItem1 = [[[NSCustomTouchBarItem alloc] initWithIdentifier:identifier] autorelease];
        self.touchBarButton1 = [[NSButton buttonWithTitle:title.toNSString() target:self
                                          action:@selector(button1Clicked)] autorelease];
        self.touchBarItem1.view =  self.touchBarButton1;
         return self.touchBarItem1;
    } else if ([identifier isEqualToString:Button2Identifier]) {
        QString title = "B2";
        self.touchBarItem2 = [[[NSCustomTouchBarItem alloc] initWithIdentifier:identifier] autorelease];
        self.touchBarButton2 = [[NSButton buttonWithTitle:title.toNSString() target:self
                                          action:@selector(button2Clicked)] autorelease];
        self.touchBarItem2.view =  self.touchBarButton2;
        return self.touchBarItem2;
    }
   return nil;
}

- (void)installAsDelegateForWindow:(NSWindow *)window
{
    _qtDelegate = window.delegate; // Save current delegate for forwarding
    window.delegate = self;
}

- (void)installAsDelegateForApplication:(NSApplication *)application
{
    _qtDelegate = application.delegate; // Save current delegate for forwarding
    application.delegate = self;
}

- (BOOL)respondsToSelector:(SEL)aSelector
{
    // We want to forward to the qt delegate. Respond to selectors it
    // responds to in addition to selectors this instance resonds to.
    return [_qtDelegate respondsToSelector:aSelector] || [super respondsToSelector:aSelector];
}

- (void)forwardInvocation:(NSInvocation *)anInvocation
{
    // Forward to the existing delegate. This function is only called for selectors
    // this instance does not responds to, which means that the qt delegate
    // must respond to it (due to the respondsToSelector implementation above).
    [anInvocation invokeWithTarget:_qtDelegate];
}

- (void)button1Clicked
{
    qDebug() << "button1Clicked";
}

- (void)button2Clicked
{
    qDebug() << "button2Clicked";
}

@end